# Script to plot residuals of the best multiple regressions for each response variable
library(nlme)
library(dplyr)
library(ggplot2)
library(patchwork)

rm(list = ls())

# 1. Load data and fit models ------------------------------------------------
load("Data/SEM_data_all_GSA.RData")

# Scale continuous variables and convert GSA to factor
SEM_data <- SEM_data %>%
  mutate(log_Length_max_cwm = log(Length_max_cwm)) %>%
  mutate(across(c(log_biomass, fric, feve, tric, teve, Length_max_cwm, log_Length_max_cwm,
                   chla_mean, depth, FPI_OTB, FPI_SSF, gravity,
                   sbt_mean, so_mean, so_sd, sst_mean, sst_sd),
                ~ as.numeric(scale(.)))) %>%
  mutate(GSA = as.factor(GSA))

# Create dummy variables for GSA (reference level is "Adriatic Sea")
SEM_data_complete <- SEM_data %>%
  # 1. Create the dummy 0/1 columns
  fastDummies::dummy_cols("GSA", remove_first_dummy = TRUE, remove_selected_columns = FALSE)

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
  random = ~1 | year, correlation = corExp(value = c(17.3, 0.2), form = ~longitude + latitude | year/GSA, nugget = TRUE), data = SEM_data_complete, method = "REML")


# 2. Plot residuals -----------------------------------------------------------
plot_diagnostics <- function(model, model_name, letter, bins = 50){

  res <- resid(model, type = "normalized")
  df <- data.frame(residuals = res, fitted = fitted(model))

  # Histogram
  p1 <- ggplot(df, aes(residuals)) +
    geom_histogram(aes(y = after_stat(density)), bins = bins, fill = "grey90", colour = "black") +
    geom_density(colour = "blue", linewidth = 1) +
    labs(title = paste0(letter, ". ", model_name), x = "Normalized residuals", y = "Density") +
    theme_classic()

  # QQ plot
  p2 <- ggplot(df, aes(sample = residuals)) +
    stat_qq(size = 1.2, alpha = 0.6) +
    stat_qq_line(colour = "blue", linewidth = 1) +
    labs(x = "Theoretical quantiles", y = "Sample quantiles") +
    theme_classic()

  # Residuals vs fitted values
  p3 <- ggplot(df, aes(x = fitted, y = residuals)) +
    geom_point(alpha = 0.5, size = 1.5) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_smooth(method = "loess", se = TRUE, colour = "blue", linewidth = 1) +
    labs(x = "Fitted values", y = "Normalized residuals") +
    theme_classic()
  
  p1 + p2 + p3 + plot_layout(ncol = 3)
}

diag_A <- plot_diagnostics(m.tric, "Taxonomic richness", "A")
diag_B <- plot_diagnostics(m.teve, "Taxonomic evenness", "B")
diag_C <- plot_diagnostics(m.fric, "Functional richness", "C")
diag_D <- plot_diagnostics(m.feve, "Functional evenness", "D")
diag_E <- plot_diagnostics(m.cwml, "Log Community-Weighted Mean maximum body length", "E")
diag_F <- plot_diagnostics(m.biomass, "Log biomass", "F")

fig_residuals <- (diag_A) /
(diag_B) /
(diag_C) /
(diag_D) /
(diag_E) /
(diag_F)

# Save figure
ggsave("Figures/Figure_S3.png", plot = fig_residuals, width = 15, height = 15, dpi = 300)
