# Script to fit SEM model for all GSAs

# Load packages
library(nlme)
library(piecewiseSEM)
library(dplyr)
library(sf)

rm(list = ls())

# 1. Load data and fit models ------------------------------------------------
load("Data/SEM_data_all_GSA.RData")

# load("Sub_model_selection/dredged_tric.RData")
# load("Sub_model_selection/dredged_teve.RData")
# load("Sub_model_selection/dredged_fric.RData")
# load("Sub_model_selection/dredged_feve.RData")
# load("Sub_model_selection/dredged_log_cwm_length.RData")
# load("Sub_model_selection/dredged_biomass.RData")

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

# Results of the best multiple regressions were adjusted by the restricted maximum likelihood method (REML) and the selected explanatory variables were retained for use in SEM
SEM_model <- psem(m.tric, m.teve, m.fric, m.feve, m.cwml, m.biomass, 
  teve %~~% tric, # added
  feve %~~% fric, # added
  feve %~~% tric, # added
  #fric %~~% teve,
  log_Length_max_cwm %~~% tric, # added
  log_Length_max_cwm %~~% fric, # added
  log_Length_max_cwm %~~% feve, # added
  log_Length_max_cwm %~~% teve, # added
  data = SEM_data_complete)

# Compute the summary of the SEM model
SEM_results <- summary(SEM_model)

rsquared(SEM_model)
SEM_results$AIC
SEM_results$dTable
SEM_results$Cstat

save(SEM_model, file = "Outputs/SEM_model_ALL_GSA.RData")
save(SEM_results, file = "Outputs/SEM_results_ALL_GSA.RData")
