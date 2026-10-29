# ==============================================================================
# Sub-model selection: functional evenness (feve)
#   Selects the predictors and the spatial correlation structure
#   component: global lme model, spatial autocorrelation (SAC) structure chosen on the full
#   model, then predictor selection by MuMIn::dredge().
#
# WARNING: long-running script (several hours)
# ==============================================================================

# ==============================================================================
# 1. Load libraries and prepare the environment
# ==============================================================================

# Load libraries
source(here::here("R/0-Load libraries.R"))

rm(list = ls())

# ==============================================================================
# 2. Data preparation, SAC structure and model selection
# ==============================================================================

options(na.action = "na.fail")   # required by dredge()

# Load the data (SEM_data)
load(here::here("Data/SEM_data_all_GSA.RData"))

# Scale continuous variables and convert GSA to factor
SEM_data <- SEM_data %>%
# Add log transformed CWM_length to the dataset
  mutate(log_Length_max_cwm = log(Length_max_cwm)) %>%
  mutate(across(c(log_biomass, fric, feve, tric, teve, Length_max_cwm, log_Length_max_cwm,
                   chla_mean, depth, FPI_OTB, FPI_SSF, gravity,
                   sbt_mean, so_mean, so_sd, sst_mean, sst_sd),
                ~ as.numeric(scale(.)))) %>%
  mutate(GSA = as.factor(GSA))

## -----------------------------------------------------------------------------
## STEP 0 - Project lon/lat to a planar CRS
## -----------------------------------------------------------------------------
sf_points <- st_as_sf(SEM_data, coords = c("longitude", "latitude"), crs = 4326, remove = FALSE)
sf_proj   <- st_transform(sf_points, crs = 3035)
proj_coords <- st_coordinates(sf_proj)

## Overwrite longitude/latitude with projected coordinates in KM
SEM_data$longitude <- proj_coords[, "X"] / 1000
SEM_data$latitude <- proj_coords[, "Y"] / 1000

cat("Coordinates projected to EPSG:3035 (ETRS89 / LAEA Europe), units = km\n")
cat("longitude (Easting)  range:", round(range(SEM_data$longitude), 1), "km\n")
cat("latitude (Northing) range:", round(range(SEM_data$latitude), 1), "km\n")

## -----------------------------------------------------------------------------
## STEP 1 - Sanity checks: collinearity among candidate (environmental) predictors
## -----------------------------------------------------------------------------

predictors <- c("chla_mean", "depth", "FPI_OTB", "FPI_SSF", "gravity",
                "sbt_mean", "so_mean", "so_sd", "sst_mean", "sst_sd",
                "teve")

cor_mat <- cor(SEM_data[, predictors], use = "pairwise.complete.obs")
print(round(cor_mat, 2))

# Flag any |r| > 0.7 for review before dredging
high_cor <- which(abs(cor_mat) > 0.7 & abs(cor_mat) < 1, arr.ind = TRUE)
if (nrow(high_cor) > 0) {
  cat("\nHighly correlated predictor pairs (|r| > 0.7):\n")
  print(unique(t(apply(high_cor, 1, function(i) sort(rownames(cor_mat)[i])))))
} else {
  cat("\nNo predictor pairs exceed |r| = 0.7\n")
}

# VIF check on environmental predictors only (plain additive check -- see
# the GVIF note below for the full model including interactions)
vif_check <- lm(as.formula(paste("feve ~", paste(predictors, collapse = " + "))),
                 data = SEM_data)
print(car::vif(vif_check))

## -----------------------------------------------------------------------------
## STEP 2 - Fit the FULL model (all environmental predictors + fixed trait x
## GSA interaction terms), ML, no spatial correlation yet
## -----------------------------------------------------------------------------

m.feve.0 <- lme(
  feve ~ chla_mean + depth + FPI_OTB + FPI_SSF + gravity + sbt_mean +
                so_mean + so_sd + sst_mean + sst_sd +
                teve * GSA,
  random = ~1 | year, data = SEM_data, method = "ML", na.action = na.fail)

summary(m.feve.0)

## -----------------------------------------------------------------------------
## STEP 3 - Estimate a starting range value for Step 4, using year x GSA cells
## -----------------------------------------------------------------------------

SEM_data$resid_feve_0 <- residuals(m.feve.0, type = "normalized")

## Build the set of year x GSA cells actually present in the data
year_gsa_cells <- unique(SEM_data[, c("year", "GSA")])

get_cell_variogram <- function(data, year_val, gsa_val, resid_col,
                                coord_x = "longitude", coord_y = "latitude") {
  sub <- data[data$year == year_val & data$GSA == gsa_val, ]
  if (nrow(sub) < 10) return(NULL)  # too few points in this cell for a reliable fit
  v <- tryCatch(
    variogram(as.formula(paste(resid_col, "~ 1")),
              locations = as.formula(paste("~", coord_x, "+", coord_y)),
              data = sub),
    error = function(e) NULL)
  if (is.null(v)) return(NULL)

  fit_converged <- TRUE
  fit <- withCallingHandlers(
    tryCatch(
      fit.variogram(v, model = vgm(psill = var(sub[[resid_col]], na.rm = TRUE),
                                    model = "Exp",
                                    range = diff(range(sub[[coord_x]], na.rm = TRUE)) / 4,
                                    nugget = 0)),
      error = function(e) NULL),
    warning = function(w) {
      fit_converged <<- FALSE
      invokeRestart("muffleWarning")
    }
  )

  list(year = year_val, gsa = gsa_val, n = nrow(sub), vgm_obj = v, fit = fit,
       converged = fit_converged)
}

cat("\nFitting variograms across", nrow(year_gsa_cells), "year x GSA cells...\n")
cell_variograms <- pbmapply(
  get_cell_variogram,
  year_val = year_gsa_cells$year,
  gsa_val  = year_gsa_cells$GSA,
  MoreArgs = list(data = SEM_data, resid_col = "resid_feve_0"),
  SIMPLIFY = FALSE
)
cell_variograms <- cell_variograms[!sapply(cell_variograms, is.null)]

n_unconverged <- sum(!sapply(cell_variograms, function(cv) cv$converged))
cat("\nUsable year x GSA cells for range estimation:", length(cell_variograms),
    "out of", nrow(year_gsa_cells), "total combinations\n")
cat(n_unconverged, "of those cells had a non-converged or singular variogram fit",
    "(excluded from the range estimate below, kept only as diagnostic info)\n")

range_estimates <- sapply(cell_variograms, function(cv) {
  if (is.null(cv$fit) || !cv$converged) return(NA)
  rng <- cv$fit$range[cv$fit$model == "Exp"]
  if (length(rng) == 1) rng else NA
})
names(range_estimates) <- sapply(cell_variograms, function(cv) paste(cv$year, cv$gsa, sep = "_"))

cat("\nPer year x GSA cell range estimates (converged fits only, for starting values):\n")
print(round(range_estimates, 2))
cat("\nMedian range estimate across converged cells:",
    round(median(range_estimates, na.rm = TRUE), 2), "km",
    "(from", sum(!is.na(range_estimates)), "converged cells)\n")

if (sum(!is.na(range_estimates)) < 5) {
  cat("\nWARNING: fewer than 5 cells produced a converged range estimate.",
      "Treat the median above with caution, or consider falling back to",
      "nlme's default correlation-structure initialization in Step 4.\n")
}

cat("\nMedian range estimate by GSA, converged cells only (context only):\n")
range_by_gsa <- tapply(range_estimates,
                        sapply(cell_variograms, function(cv) cv$gsa),
                        median, na.rm = TRUE)
print(round(range_by_gsa, 2))

## -----------------------------------------------------------------------------
## STEP 4 - Choose SAC structure on the FULL model (before any dredging)
## -----------------------------------------------------------------------------

start_range <- if (exists("range_estimates")) median(range_estimates, na.rm = TRUE) else NA
if (is.na(start_range)) start_range <- NULL
cat("\nUsing starting range value for correlation structures:",
    ifelse(is.null(start_range), "nlme default", round(start_range, 2)), "\n")

init_val <- if (!is.null(start_range)) c(start_range, 0.2) else NULL

spatial_form <- ~longitude + latitude | year/GSA

fit_one <- function(base_model, cor_fun_name, form, use_value, init_value, label) {
  tryCatch({
    if (use_value) {
      cor_call <- bquote(
        .(as.name(cor_fun_name))(value = .(init_value), form = .(form), nugget = TRUE)
      )
    } else {
      cor_call <- bquote(
        .(as.name(cor_fun_name))(form = .(form), nugget = TRUE)
      )
    }
    m <- eval(bquote(update(.(base_model), correlation = .(cor_call))))
    m
  }, error = function(e) {
    cat(label, "failed:", conditionMessage(e), "\n")
    NULL
  })
}

cor_families <- list(exp = "corExp", gaus = "corGaus", spher = "corSpher", ratio = "corRatio", lin = "corLin")
takes_value  <- list(exp = TRUE,   gaus = TRUE,    spher = TRUE,     ratio = TRUE,     lin = FALSE)

candidate_models <- list(none = m.feve.0)

cat("\nFitting", length(cor_families), "correlation structures on the full model",
    "(pre-dredge)...\n")
pb <- txtProgressBar(min = 0, max = length(cor_families), style = 3)
for (i in seq_along(cor_families)) {
  fam_name <- names(cor_families)[i]
  candidate_models[[fam_name]] <- fit_one(
    base_model   = m.feve.0,
    cor_fun_name = cor_families[[fam_name]],
    form         = spatial_form,
    use_value    = takes_value[[fam_name]],
    init_value   = init_val,
    label        = fam_name
  )
  setTxtProgressBar(pb, i)
}
close(pb)

candidate_models <- candidate_models[!sapply(candidate_models, is.null)]

aic_table <- data.frame(structure = names(candidate_models), AIC = sapply(candidate_models, AIC))
aic_table <- aic_table[order(aic_table$AIC), ]
print(aic_table)

best_struct_name <- aic_table$structure[1]
cat("\nBest structure (chosen on the full, pre-dredge model):", best_struct_name, "\n")
if (best_struct_name != "none") {
  print(summary(candidate_models[[best_struct_name]])$modelStruct$corStruct)
}

m.feve.sac <- candidate_models[[best_struct_name]]

## -----------------------------------------------------------------------------
## STEP 5 - Dredge the ENVIRONMENTAL predictors only, under the fixed SAC
## structure. Trait variables and their GSA interactions stay in every
## candidate model via dredge() fixed argument
## -----------------------------------------------------------------------------

tl <- attr(terms(m.feve.sac), "term.labels")
print(tl)  # inspect, then copy the 5 interaction terms exactly as shown

fixed_terms <- c("teve", "GSA", tl[grepl("GSA", tl) & grepl("teve", tl)])

fixed_terms_formula <- as.formula(paste("~", paste(fixed_terms, collapse = " + ")))

# Warning: this step can take several hours to complete !!!!
dredged.feve <- dredge(m.feve.sac, fixed = fixed_terms_formula, trace = 2)

# dir.create(here::here("Outputs/Sub_model_selection"), showWarnings = FALSE, recursive = TRUE)
# save(dredged.feve, file = here::here("Outputs/Sub_model_selection/dredged_feve.RData"))

print(head(model.sel(dredged.feve), 10))
print(sw(dredged.feve))


top_set <- subset(dredged.feve, delta < 2)
cat("\nNumber of models within 2 AICc units of the best:", nrow(top_set), "\n")
 
m.feve.final.ml <- get.models(dredged.feve, subset = 1)[[1]]
summary(m.feve.final.ml)

## -----------------------------------------------------------------------------
## STEP 6 - Finalize with REML + residual diagnostics
## -----------------------------------------------------------------------------

m.feve.final <- update(m.feve.final.ml, method = "REML")
summary(m.feve.final)

par(mfrow = c(1, 3))
hist(residuals(m.feve.final, type = "normalized"), main = "Histogram of residuals",
     xlab = "Normalized residuals")
plot(fitted(m.feve.final), residuals(m.feve.final, type = "normalized"),
     xlab = "Fitted", ylab = "Normalized residuals", main = "Residuals vs Fitted")
abline(h = 0, lty = 2, col = "red")
qqnorm(residuals(m.feve.final, type = "normalized"), main = "QQ plot")
qqline(residuals(m.feve.final, type = "normalized"), col = "red")
par(mfrow = c(1, 1))

## Re-check spatial autocorrelation on FINAL residuals, per year x GSA cell
SEM_data$resid_feve_final <- residuals(m.feve.final, type = "normalized")

check_cell_autocorr <- function(data, year_val, gsa_val, resid_col) {
  sub <- data[data$year == year_val & data$GSA == gsa_val, ]
  if (nrow(sub) < 10) return(invisible(NULL))
  dist_mat <- as.matrix(dist(cbind(sub$longitude, sub$latitude)))
  inv_dist <- 1 / dist_mat
  diag(inv_dist) <- 0
  inv_dist[is.infinite(inv_dist)] <- 0
  mi <- tryCatch(Moran.I(sub[[resid_col]], inv_dist), error = function(e) NULL)
  if (is.null(mi)) return(invisible(NULL))
  cat(sprintf("Year %s | GSA %s - Moran's I = %.3f, p = %.4f (n = %d)\n",
              year_val, gsa_val, mi$observed, mi$p.value, nrow(sub)))
  invisible(mi)
}

cat("\n--- Spatial autocorrelation check per year x GSA cell (FINAL residuals) ---\n")
invisible(mapply(check_cell_autocorr,
                  year_val = year_gsa_cells$year,
                  gsa_val  = year_gsa_cells$GSA,
                  MoreArgs = list(data = SEM_data, resid_col = "resid_feve_final")))
