library(nlme)
library(piecewiseSEM)
library(marginaleffects)

rm(list = ls())

# 1. Load data 
load("Data/SEM_data_all_GSA.RData")

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

# Define the biomass model
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


### Refit log biomass sub-model with taxonomic and functional metrics separately
# drops every term containing a metric (main effect + GSA interaction), factor or dummy coding
drop_metrics <- function(model, metrics) {
  lab <- attr(terms(formula(model)), "term.labels")
  hit <- sapply(strsplit(lab, ":", fixed = TRUE), function(p) any(p %in% metrics))
  update(model, fixed. = as.formula(paste(". ~ . -", paste(lab[hit], collapse = " - "))))
}

# Define the models to compare
m_ML <- update(m.biomass, method = "ML")   # ML requires to compare different fixed effects

mods <- list(
  main     = m_ML,                                   # tric + teve + fric + feve + CWM
  no_fric  = drop_metrics(m_ML, "fric"),             # taxonomic richness only
  no_tric  = drop_metrics(m_ML, "tric"),             # functional richness only
  tax_only = drop_metrics(m_ML, c("fric", "feve")),  # taxonomic metrics only
  fun_only = drop_metrics(m_ML, c("tric", "teve"))   # functional metrics only
)

cmp <- data.frame(model = names(mods),
  df  = sapply(mods, function(m) attr(logLik(m), "df")),
  AIC = sapply(mods, AIC), BIC = sapply(mods, BIC),
  R2m = sapply(mods, function(m) rsquared(m)$Marginal),
  R2c = sapply(mods, function(m) rsquared(m)$Conditional))

cmp$dAIC <- cmp$AIC - min(cmp$AIC)
print(cmp)

anova(mods$no_fric, mods$main)   # does functional richness add to taxonomic richness?
anova(mods$no_tric, mods$main)   # does taxonomic richness add to functional richness?
anova(mods$tax_only, mods$main)  # do functional metrics add to taxonomic metrics?
anova(mods$fun_only, mods$main)  # do taxonomic metrics add to functional metrics?

### Slopes by GSA across specifications and collinearity
metrics <- c("tric","teve","fric","feve","log_Length_max_cwm")

slopes_by_model <- dplyr::bind_rows(lapply(names(mods), function(n) {
  lab <- attr(terms(formula(mods[[n]])), "term.labels")
  avg_slopes(mods[[n]], variables = intersect(metrics, lab), by = "GSA", newdata = SEM_data_complete) %>%
    as.data.frame() %>% dplyr::mutate(model = n)
}))

slopes_by_model %>% dplyr::filter(term == "tric") %>%
  dplyr::select(model, GSA, estimate, std.error, p.value)

slopes_by_model %>% dplyr::filter(term == "fric") %>%
  dplyr::select(model, GSA, estimate, std.error, p.value)

by(SEM_data_complete[, c("tric","fric")], SEM_data_complete$GSA, function(d) cor(d$tric, d$fric))
car::vif(lm(log_biomass ~ tric + fric + teve + feve + log_Length_max_cwm + GSA, data = SEM_data_complete))

### Predictive performance with spatial block cross-validation
dat <- as.data.frame(SEM_data_complete); set.seed(1)
dat$haul_id <- seq_len(nrow(dat))

# spatial blocks within each GSA, assigned to 5 folds stratified by GSA
dat$block <- NA_character_
for (g in levels(dat$GSA)) {
  i <- which(dat$GSA == g)
  dat$block[i] <- paste(g, kmeans(dat[i, c("longitude","latitude")], centers = 10)$cluster, sep = "_")
}
blk <- unique(dat[, c("block","GSA")])
blk$fold <- ave(seq_len(nrow(blk)), blk$GSA, FUN = function(i) sample(rep_len(1:5, length(i))))
dat$fold <- blk$fold[match(dat$block, blk$block)]

# CV for one model: returns out-of-fold predictions + one row of scores per fold
cv_perf <- function(model, data, fold_col = "fold", response = "log_biomass") {
  pred <- rep(NA_real_, nrow(data)); fs <- list()
  for (k in sort(unique(data[[fold_col]]))) {
    idx <- data[[fold_col]] == k
    tr <- data[!idx, ]; te <- data[idx, ]
    fit <- tryCatch(update(model, data = tr),
             error = function(e) tryCatch(update(model, data = tr,
               control = lmeControl(opt = "optim", maxIter = 200, msMaxIter = 200)),
               error = function(e2) NULL))
    if (is.null(fit)) {
      fs[[as.character(k)]] <- data.frame(fold = k, n_test = nrow(te), RMSE = NA_real_, R2 = NA_real_, converged = FALSE)
      next
    }
    p <- predict(fit, newdata = te, level = 0); y <- te[[response]]
    pred[idx] <- p
    fs[[as.character(k)]] <- data.frame(fold = k, n_test = nrow(te), RMSE = sqrt(mean((y - p)^2)),
                                        R2 = 1 - sum((y - p)^2) / sum((y - mean(y))^2), converged = TRUE)
  }
  list(pred = pred, folds = do.call(rbind, fs))
}
cv_res <- lapply(mods, cv_perf, data = dat)      # slow step: 5 models x 5 folds

score <- function(y, p) { ok <- !is.na(p)
  c(RMSE = sqrt(mean((y[ok] - p[ok])^2)),
    R2   = 1 - sum((y[ok] - p[ok])^2) / sum((y[ok] - mean(y[ok]))^2)) }

# (a) main table: one row per model
cv_summary <- do.call(rbind, lapply(names(cv_res), function(n) {
  s <- score(dat$log_biomass, cv_res[[n]]$pred); f <- cv_res[[n]]$folds
  data.frame(model = n, RMSE_cv = s["RMSE"], R2_cv = s["R2"],
             RMSE_fold_mean = mean(f$RMSE, na.rm = TRUE), RMSE_fold_sd = sd(f$RMSE, na.rm = TRUE),
             R2_fold_mean   = mean(f$R2,   na.rm = TRUE), R2_fold_sd   = sd(f$R2,   na.rm = TRUE),
             folds_failed = sum(!f$converged), row.names = NULL)
}))

# (b) per-fold scores and paired differences against the main model
cv_folds <- do.call(rbind, lapply(names(cv_res), function(n) cbind(model = n, cv_res[[n]]$folds)))
ref <- cv_folds[cv_folds$model == "main", c("fold", "RMSE", "R2")]
names(ref) <- c("fold", "RMSE_main", "R2_main")
cv_delta <- merge(cv_folds, ref, by = "fold")
cv_delta$dRMSE <- cv_delta$RMSE - cv_delta$RMSE_main   # > 0: worse than main
cv_delta$dR2   <- cv_delta$R2   - cv_delta$R2_main     # < 0: worse than main
cv_delta_summary <- do.call(rbind, lapply(split(cv_delta, cv_delta$model), function(d)
  data.frame(model = d$model[1], mean_dRMSE = mean(d$dRMSE, na.rm = TRUE),
             mean_dR2 = mean(d$dR2, na.rm = TRUE),
             folds_better_than_main = sum(d$dRMSE < 0, na.rm = TRUE), n_folds = sum(!is.na(d$dRMSE)))))

# (c) performance by GSA
cv_by_gsa <- do.call(rbind, lapply(names(cv_res), function(n)
  do.call(rbind, lapply(levels(dat$GSA), function(g) {
    i <- dat$GSA == g; s <- score(dat$log_biomass[i], cv_res[[n]]$pred[i])
    data.frame(model = n, GSA = g, RMSE_cv = s["RMSE"], R2_cv = s["R2"], row.names = NULL) }))))

# (d) out-of-fold predictions, one row per haul (for observed-vs-predicted plots)
cv_oof <- data.frame(haul_id = dat$haul_id, GSA = dat$GSA, year = dat$year, block = dat$block,
                     fold = dat$fold, observed = dat$log_biomass, sapply(cv_res, function(x) x$pred))
names(cv_oof)[-(1:6)] <- paste0("pred_", names(cv_res))

# ---- save ----
out_dir <- "Sensitivity_analyses"; dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
saveRDS(list(cv_res = cv_res, fold_assignment = dat[, c("haul_id","GSA","block","fold")],
             summary = cv_summary, folds = cv_folds, delta = cv_delta,
             by_gsa = cv_by_gsa, oof = cv_oof), file.path(out_dir, "cv_results.rds"))
write.csv(cv_summary,       file.path(out_dir, "cv_summary.csv"),       row.names = FALSE)
write.csv(cv_folds,         file.path(out_dir, "cv_per_fold.csv"),      row.names = FALSE)
write.csv(cv_delta_summary, file.path(out_dir, "cv_delta_vs_main.csv"), row.names = FALSE)
write.csv(cv_by_gsa,        file.path(out_dir, "cv_by_gsa.csv"),        row.names = FALSE)
write.csv(cv_oof,           file.path(out_dir, "cv_out_of_fold_predictions.csv"), row.names = FALSE)

