# ==============================================================================
#Residual spatial autocorrelation check of the final sub-models
#
#   Refits the six final sub-models (with and without their spatial correlation structure)
#   and tests residual spatial autocorrelation (Moran's I, inverse-distance weights) for each
#   year x GSA cell, with Benjamini-Hochberg correction.
#
# ==============================================================================

# ==============================================================================
# 1. Load libraries and prepare the environment
# ==============================================================================

# Load libraries
source(here::here("R/0-Load libraries.R"))

rm(list = ls())

# 1. Load data and fit models ------------------------------------------------

load(here::here("Data/SEM_data_all_GSA.RData"))

# load(here::here("Outputs/Sub_model_selection/dredged_tric.RData"))
# load(here::here("Outputs/Sub_model_selection/dredged_teve.RData"))
# load(here::here("Outputs/Sub_model_selection/dredged_fric.RData"))
# load(here::here("Outputs/Sub_model_selection/dredged_feve.RData"))
# load(here::here("Outputs/Sub_model_selection/dredged_log_cwm_length.RData"))
# load(here::here("Outputs/Sub_model_selection/dredged_biomass.RData"))

# Prepare the data
SEM_data <- SEM_data %>%
  mutate(
    log_Length_max_cwm = log(Length_max_cwm),
    across(c(log_biomass, fric, feve, tric, teve, Length_max_cwm, log_Length_max_cwm,
             chla_mean, depth, FPI_OTB, FPI_SSF, gravity, sbt_mean, so_mean, so_sd,
             sst_mean, sst_sd), ~ as.numeric(scale(.))),
    GSA = as.factor(GSA)
  )

# Project coordinates (km)
proj_coords <- SEM_data %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326, remove = FALSE) %>%
  st_transform(3035) %>%
  st_coordinates()

SEM_data[c("longitude", "latitude")] <- proj_coords / 1000

# Create dummy variables for GSA (reference level is "Adriatic Sea")
SEM_data_complete <- SEM_data %>%
  # 1. Create the dummy 0/1 columns
  fastDummies::dummy_cols("GSA", remove_first_dummy = TRUE, remove_selected_columns = FALSE)

# Define sub-models based on previous analyses
m.tric <- nlme::lme(tric ~ chla_mean + depth + FPI_OTB + FPI_SSF + gravity + sbt_mean + so_mean + so_sd + sst_sd +
  # GSA main effects
         GSA_Aegean_Sea + GSA_Catalan_Sea + GSA_Gulf_of_Lion,
  random = ~1 | year, correlation = corExp(value = c(37.1, 0.2), form = ~longitude + latitude | year/GSA, nugget = TRUE),
  data = SEM_data_complete, method = "REML", na.action = na.fail)

m.teve <- nlme::lme(teve ~ chla_mean + sbt_mean + so_mean + so_sd + sst_mean + sst_sd +
  # GSA main effects
         GSA_Aegean_Sea + GSA_Catalan_Sea + GSA_Gulf_of_Lion,
  random = ~1 | year, correlation = corRatio(value = c(19.1, 0.2), form = ~longitude + latitude | year/GSA, nugget = TRUE),
  data = SEM_data_complete, method = "REML", na.action = na.fail)

m.fric <- nlme::lme(fric ~ chla_mean + FPI_OTB + FPI_SSF + gravity + sbt_mean + so_sd + sst_mean + sst_sd +
         # GSA main effects
         GSA_Aegean_Sea + GSA_Catalan_Sea + GSA_Gulf_of_Lion +
         # tric main effect
         tric +
         # tric:GSA interactions
         tric:GSA_Aegean_Sea + tric:GSA_Catalan_Sea + tric:GSA_Gulf_of_Lion,
  random = ~1 | year, correlation = corRatio(value = c(16.7, 0.2), form = ~longitude + latitude | year/GSA, nugget = TRUE),
  data = SEM_data_complete, method = "REML", na.action = na.fail)

m.feve <- nlme::lme(
  feve ~ chla_mean + depth + gravity + sbt_mean + sst_mean +
         # GSA main effects
         GSA_Aegean_Sea + GSA_Catalan_Sea + GSA_Gulf_of_Lion +
         # teve main effect
         teve +
         # teve:GSA interactions
         teve:GSA_Aegean_Sea + teve:GSA_Catalan_Sea + teve:GSA_Gulf_of_Lion,
  random = ~1 | year, correlation = corExp(value = c(13.5, 0.2), form = ~longitude + latitude | year/GSA, nugget = TRUE),
  data = SEM_data_complete, method = "REML", na.action = na.fail)

m.cwml <- nlme::lme(log_Length_max_cwm ~  chla_mean +depth + FPI_SSF + gravity + sbt_mean + so_mean + sst_mean + sst_sd + 
  # GSA main effects
         GSA_Aegean_Sea + GSA_Catalan_Sea + GSA_Gulf_of_Lion,
  random = ~1 | year, correlation = corExp(value = c(16.2, 0.2), form = ~longitude + latitude | year/GSA, nugget = TRUE),
  data = SEM_data_complete, method = "REML", na.action = na.fail)

m.biomass <- nlme::lme(log_biomass ~
         # environment main effects
        chla_mean + depth + FPI_OTB + sbt_mean + so_mean + so_sd + sst_mean +
         # GSA main effects
         GSA_Aegean_Sea + GSA_Catalan_Sea + GSA_Gulf_of_Lion +
         # trait/diversity main effects
         fric + feve + log_Length_max_cwm + teve + tric +
         # fric:GSA interactions
         fric:GSA_Aegean_Sea + fric:GSA_Catalan_Sea + fric:GSA_Gulf_of_Lion +
         # feve:GSA interactions
         feve:GSA_Aegean_Sea + feve:GSA_Catalan_Sea + feve:GSA_Gulf_of_Lion +
         # Length_max_cwm:GSA interactions
         log_Length_max_cwm:GSA_Aegean_Sea + log_Length_max_cwm:GSA_Catalan_Sea + log_Length_max_cwm:GSA_Gulf_of_Lion +
         # teve:GSA interactions
         teve:GSA_Aegean_Sea + teve:GSA_Catalan_Sea + teve:GSA_Gulf_of_Lion +
         # tric:GSA interactions
         tric:GSA_Aegean_Sea + tric:GSA_Catalan_Sea + tric:GSA_Gulf_of_Lion,
  random = ~1 | year, correlation = corExp(value = c(17.3, 0.2), form = ~longitude + latitude | year/GSA, nugget = TRUE),
  data = SEM_data_complete, method = "REML", na.action = na.fail)


## -----------------------------------------------------------------------------
# ==============================================================================
# 2b. Moran's I of the residuals
# ==============================================================================

## 1. Final sub-models (with spatial correlation structure)
## -----------------------------------------------------------------------------
models <- list(
  tric    = m.tric,
  teve    = m.teve,
  fric    = m.fric,
  feve    = m.feve,
  cwml    = m.cwml,
  biomass = m.biomass
)

## -----------------------------------------------------------------------------
## 2. Moran's I for one year x GSA cell (inverse-distance weights)
## -----------------------------------------------------------------------------
moran_cell <- function(dat, resid_col, y, g, min_n = 10) {
  sub <- dat[dat$year == y & dat$GSA == g, ]
  if (nrow(sub) < min_n) return(NULL)
  d <- as.matrix(dist(cbind(sub$longitude, sub$latitude)))   # km (projected)
  w <- 1 / d
  diag(w) <- 0
  w[is.infinite(w)] <- 0                                     # duplicated coordinates
  mi <- ape::Moran.I(sub[[resid_col]], w)
  data.frame(year = as.character(y), GSA = as.character(g), n = nrow(sub),
             I = mi$observed, E = mi$expected, p = mi$p.value)
}

## -----------------------------------------------------------------------------
## 3. Moran's I over all groups (gsa x year) for one fitted model
## -----------------------------------------------------------------------------
moran_model <- function(m, label, version) {
  dat <- nlme::getData(m)                                    
  dat$res <- residuals(m, type = "normalized")
  dat$GSA <- as.character(dat$GSA)
  cells <- unique(dat[, c("year", "GSA")])
  out <- do.call(rbind, Map(function(y, g) moran_cell(dat, "res", y, g),
                            cells$year, cells$GSA))
  out$model   <- label
  out$version <- version
  out
}

## -----------------------------------------------------------------------------
## 4. Final models (with correlation) and the same models without it
## -----------------------------------------------------------------------------
res_after <- Map(function(m, nm) moran_model(m, nm, "with spatial structure"),
                 models, names(models))

models_nocor <- lapply(models, function(m) update(m, correlation = NULL))
res_before <- Map(function(m, nm) moran_model(m, nm, "no spatial structure"),
                  models_nocor, names(models_nocor))

moran_cells <- do.call(rbind, c(unname(res_before), unname(res_after)))
rownames(moran_cells) <- NULL

## Check the structure before grouping
str(moran_cells)
table(moran_cells$model, moran_cells$version)

moran_cells <- moran_cells %>%
  dplyr::group_by(model, version) %>%
  dplyr::mutate(p_BH = p.adjust(p, method = "BH")) %>%   # BH across groups within each model x version
  dplyr::ungroup() %>%
  dplyr::select(model, version, GSA, year, n, I, E, p, p_BH) %>%
  dplyr::arrange(model, version, GSA, year)

## -----------------------------------------------------------------------------
## 5. Summaries
## -----------------------------------------------------------------------------
summ_fun <- function(df) {
  summarise(df,
            cells    = n(),
            sig_raw  = sum(p < 0.05),
            sig_BH   = sum(p_BH < 0.05),
            pct_sig  = round(100 * mean(p < 0.05), 0),
            median_I = round(median(I), 3),
            min_I    = round(min(I), 3),
            max_I    = round(max(I), 3),
            .groups = "drop")
}

moran_by_gsa <- moran_cells %>%
  group_by(model, version, GSA) %>% summ_fun()

moran_overall <- moran_cells %>%
  group_by(model, version) %>% summ_fun() %>%
  mutate(GSA = "All", .after = version)

moran_summary <- bind_rows(moran_by_gsa, moran_overall) %>%
  arrange(model, GSA, version)

## Wide version for the appendix table (before vs after side by side)
moran_wide <- moran_summary %>%
  mutate(version = ifelse(version == "no spatial structure", "before", "after")) %>%
  pivot_wider(id_cols = c(model, GSA), names_from = version,
              values_from = c(cells, sig_raw, sig_BH, median_I, min_I, max_I),
              names_glue = "{.value}_{version}")

## -----------------------------------------------------------------------------
## 6. Export
## -----------------------------------------------------------------------------
# write_xlsx(list(cells   = moran_cells,
          #       summary = moran_summary,
          #       wide    = moran_wide),
          #  path = here::here("Outputs/Moran_I_residuals.xlsx"))
