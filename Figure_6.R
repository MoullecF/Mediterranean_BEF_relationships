library(nlme)
library(marginaleffects)
library(dplyr)
library(ggplot2)

rm(list = ls())

# Load the data
load("Data/SEM_data_all_GSA.RData")

# Data preprocessing: scale numeric variables and set GSA as a factor with specific levels
SEM_data <- SEM_data %>%
# Add log transformed CWM_length to the dataset
  mutate(log_Length_max_cwm = log(Length_max_cwm)) %>%
  mutate(across(c(log_biomass, fric, feve, tric, teve, Length_max_cwm, log_Length_max_cwm,
                   chla_mean, depth, FPI_OTB, FPI_SSF, gravity,
                   sbt_mean, so_mean, so_sd, sst_mean, sst_sd),
                ~ as.numeric(scale(.)))) %>%
  mutate(GSA = as.factor(GSA)) %>%
  mutate(GSA = factor(GSA, levels = c("Adriatic_Sea", "Aegean_Sea", "Catalan_Sea", "Gulf_of_Lion")))

# Project coordinates (km)
proj_coords <- SEM_data %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326, remove = FALSE) %>%
  st_transform(3035) %>%
  st_coordinates()

SEM_data[c("longitude", "latitude")] <- proj_coords / 1000

# Fit the linear mixed-effects model with interactions between biodiversity metrics and GSA
m.biomass <- nlme::lme(
  log_biomass ~ chla_mean + depth + FPI_OTB + sbt_mean + so_mean + so_sd + sst_mean +
                (fric + log_Length_max_cwm + teve + tric + feve) * GSA,
                random = ~1 | year, correlation = corExp(value = c(17.3, 0.2), form = ~longitude + latitude | year/GSA, nugget = TRUE),
                data = SEM_data, method = "REML")

# ==============================================================================
# 1. ANOVA 
# ==============================================================================
# Patch for MuMIn::model.frame.lme to handle random effects correctly
if ("MuMIn" %in% loadedNamespaces()) {
  model.frame.lme_fixed <- function(formula, data, random = FALSE, ...) {
    x <- formula
    frm <- formula(x)
    if (isTRUE(random)) {
      for (reStruct in x$modelStruct$reStruct) {
        frm[[3L]] <- call("+", frm[[3L]], attr(reStruct, "formula")[[2L]])
      }
    }
    if (missing(data)) data <- x$data[rownames(x$fitted), ]
    do.call("model.frame", list(formula = frm, data = data, drop.unused.levels = TRUE))
  }
  assignInNamespace("model.frame.lme", model.frame.lme_fixed, "MuMIn")
}

# Perform an ANOVA to test the significance of the model terms
omnibus_test <- car::Anova(m.biomass, type = "III")
print(omnibus_test) 

# ==============================================================================
# 2. PAIRWISE COMPARISONS OF SLOPES BETWEEN GSA
# ==============================================================================
compare_slopes_me <- function(model, var, metric_name, raw_data = SEM_data) {
  
contrasts_gsa <- marginaleffects::avg_slopes(
    model, variables = var, by = "GSA", newdata = raw_data,
    hypothesis = ~ pairwise
  ) %>%
    as.data.frame() %>%
    dplyr::mutate(
      metric = metric_name,
      p.value = stats::p.adjust(p.value, method = "holm"),
      stars = dplyr::case_when(
        p.value <= 0.001 ~ "***",
        p.value <= 0.01  ~ "**",
        p.value <= 0.05  ~ "*",
        p.value <= 0.1   ~ ".",
        TRUE             ~ "ns"
      )
    )
  
  return(contrasts_gsa)
}

# Combine all biodiversity metrics
slope_comparisons_me <- dplyr::bind_rows(
  compare_slopes_me(m.biomass, "feve", "Functional evenness"),
  compare_slopes_me(m.biomass, "fric", "Functional richness"),
  compare_slopes_me(m.biomass, "tric", "Taxonomic richness"),
  compare_slopes_me(m.biomass, "teve", "Taxonomic evenness"),
  compare_slopes_me(m.biomass, "log_Length_max_cwm", "Log CWM length")
)

# ==============================================================================
# 3. CLEAN STRINGS AND PLOT MAPPING
# ==============================================================================
slope_comparisons_me <- slope_comparisons_me %>%
  dplyr::mutate(contrast_clean = hypothesis %>% 
    stringr::str_remove_all("[\\(\\)]") %>% # Retire les parenthèses
    stringr::str_replace_all("_", " ") %>%  # Remplace les underscores par des espaces
    stringr::str_replace(" - ", " vs ")    # Transforme le moins en "vs"
  ) %>%
  dplyr::mutate(
    flip = dplyr::case_when(
      stringr::str_detect(contrast_clean, "Aegean Sea vs Adriatic Sea")   ~ TRUE,
      stringr::str_detect(contrast_clean, "Catalan Sea vs Adriatic Sea")  ~ TRUE,
      stringr::str_detect(contrast_clean, "Gulf of Lion vs Adriatic Sea") ~ TRUE,
      stringr::str_detect(contrast_clean, "Catalan Sea vs Aegean Sea")   ~ TRUE,
      stringr::str_detect(contrast_clean, "Gulf of Lion vs Aegean Sea")  ~ TRUE,
      stringr::str_detect(contrast_clean, "Gulf of Lion vs Catalan Sea") ~ TRUE,
      TRUE ~ FALSE
    ),
    statistic = dplyr::if_else(flip, -statistic, statistic),
    conf.low_new  = dplyr::if_else(flip, -conf.high, conf.low),
    conf.high_new = dplyr::if_else(flip, -conf.low, conf.high),
    conf.low  = conf.low_new,
    conf.high = conf.high_new,
    estimate = dplyr::if_else(flip, -estimate, estimate)
  ) %>%
  dplyr::select(-flip, -conf.low_new, -conf.high_new) %>%
  dplyr::mutate(
    contrast_final = dplyr::case_when(
      contrast_clean == "Aegean Sea vs Adriatic Sea"   ~ "Adriatic Sea vs Aegean Sea",
      contrast_clean == "Catalan Sea vs Adriatic Sea"  ~ "Adriatic Sea vs Northern Spain",
      contrast_clean == "Gulf of Lion vs Adriatic Sea" ~ "Adriatic Sea vs Gulf of Lion",
      contrast_clean == "Catalan Sea vs Aegean Sea"   ~ "Aegean Sea vs Northern Spain",
      contrast_clean == "Gulf of Lion vs Aegean Sea"  ~ "Aegean Sea vs Gulf of Lion",
      contrast_clean == "Gulf of Lion vs Catalan Sea" ~ "Northern Spain vs Gulf of Lion",
      TRUE ~ contrast_clean
    )
  ) %>%
  # Set the order of contrasts for plotting
  dplyr::mutate(contrast_final = factor(contrast_final, levels = c(
    "Adriatic Sea vs Aegean Sea", "Adriatic Sea vs Northern Spain", "Adriatic Sea vs Gulf of Lion",
    "Aegean Sea vs Northern Spain", "Aegean Sea vs Gulf of Lion", "Northern Spain vs Gulf of Lion"
  ))) %>%
# order biodiversity metrics for plotting
  dplyr::mutate(metric = factor(metric, levels = c("Log CWM length", "Functional evenness", "Functional richness", "Taxonomic evenness", "Taxonomic richness")))

# Calculate the maximum absolute value of the slope estimates for color scaling
max_abs <- max(abs(slope_comparisons_me$estimate), na.rm = TRUE)

# Create the slope comparison plot
slope_plot <- ggplot(slope_comparisons_me, aes(x = contrast_final, y = metric, fill = estimate)) +
  geom_tile(color = "white") +
  geom_text(aes(label = stars), size = 4, fontface = "bold") +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, limits = c(-max_abs, max_abs), name = "Slope\ndifference") +
  theme_classic(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title = element_text(face = "bold", size = 11)) +
  labs(x = "GSA comparison", y = "Biodiversity metrics", title = "Pairwise differences between GSA slopes")

ggsave("Figures/Figure_6.png", plot = slope_plot, width = 6, height = 4, dpi = 300)
