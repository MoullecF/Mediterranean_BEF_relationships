# =========================
# 1. LOAD PACKAGES
# =========================
library(nlme)
library(ggeffects)
library(dplyr)
library(ggplot2)
library(patchwork)
library(sf)

rm(list = ls())

load("Data/SEM_data_all_GSA.RData")

# =========================
# 2. PREPARE DATA
# =========================
SEM_data <- SEM_data %>%
# Add log transformed CWM_length to the dataset
  mutate(log_Length_max_cwm = log(Length_max_cwm)) %>%
  mutate(across(c(log_biomass, fric, feve, tric, teve, Length_max_cwm, log_Length_max_cwm,
                   chla_mean, depth, FPI_OTB, FPI_SSF, gravity,
                   sbt_mean, so_mean, so_sd, sst_mean, sst_sd),
                ~ as.numeric(scale(.)))) %>%
  mutate(GSA = as.factor(GSA))

# Project coordinates (km)
proj_coords <- SEM_data %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326, remove = FALSE) %>%
  st_transform(3035) %>%
  st_coordinates()

SEM_data[c("longitude", "latitude")] <- proj_coords / 1000


data_long <- SEM_data %>%
  select(log_biomass, GSA, sst_mean, sst_sd, so_mean, so_sd, sbt_mean, chla_mean, depth, FPI_OTB, FPI_SSF, gravity) %>%
  tidyr::pivot_longer(
    cols = c(sst_mean, sst_sd, so_mean, so_sd, sbt_mean, chla_mean, depth, FPI_OTB, FPI_SSF, gravity),
    names_to = "metric",
    values_to = "value") %>%
  mutate(metric = recode(metric,
    sst_mean = "SST mean",
    sst_sd   = "SST variability",
    so_mean  = "Salinity mean",
    so_sd    = "Salinity variability",
    sbt_mean = "SBT mean",
    chla_mean = "Chlorophyll-a mean",
    depth = "Depth",
    FPI_OTB = "Fishing pressure (OTB)",
    FPI_SSF = "Fishing pressure (SSF)",
    gravity = "Gravity")) %>%
  mutate(metric = factor(metric, levels = c("SST mean", "SST variability", "Salinity mean", "Salinity variability", "SBT mean", "Chlorophyll-a mean", "Depth", "Fishing pressure (OTB)", "Fishing pressure (SSF)", "Gravity")))

# =========================
# 3. FIT MODEL
# =========================
m.biomass <- nlme::lme(
  log_biomass ~ chla_mean + depth + FPI_OTB + sbt_mean + so_mean + so_sd + sst_mean +
                (fric + log_Length_max_cwm + teve + tric + feve) * GSA,
                random = ~1 | year, correlation = corExp(value = c(17.3, 0.2), form = ~longitude + latitude | year/GSA, nugget = TRUE),
                data = SEM_data, method = "REML")

# =========================
# 4. GENERATE PREDICTIONS
# =========================
pred_sst <- ggpredict(m.biomass, terms = c("sst_mean [all]", "GSA")) %>%
  as.data.frame() %>%
  mutate(metric = "SST mean")

# pred_sst_sd <- ggpredict(m.biomass, terms = c("sst_sd [all]", "GSA")) %>%
#   as.data.frame() %>%
#   mutate(metric = "SST variability") 

pred_so <- ggpredict(m.biomass, terms = c("so_mean [all]", "GSA")) %>%
  as.data.frame() %>%
  mutate(metric = "Salinity mean")

pred_so_sd <- ggpredict(m.biomass, terms = c("so_sd [all]", "GSA")) %>%
  as.data.frame() %>%
  mutate(metric = "Salinity variability")

pred_chla <- ggpredict(m.biomass, terms = c("chla_mean [all]", "GSA")) %>%
  as.data.frame() %>%
  mutate(metric = "Chlorophyll-a mean")

pred_sbt <- ggpredict(m.biomass, terms = c("sbt_mean [all]", "GSA")) %>%
  as.data.frame() %>%
  mutate(metric = "SBT mean")

pred_depth <- ggpredict(m.biomass, terms = c("depth [all]", "GSA")) %>%
  as.data.frame() %>%
  mutate(metric = "Depth")

pred_fpi_otb <- ggpredict(m.biomass, terms = c("FPI_OTB [all]", "GSA")) %>%
  as.data.frame() %>%
  mutate(metric = "Fishing pressure (OTB)")

# pred_fpi_ssf <- ggpredict(m.biomass, terms = c("FPI_SSF [all]", "GSA")) %>%
#   as.data.frame() %>%
#   mutate(metric = "Fishing pressure (SSF)")

# pred_gravity <- ggpredict(m.biomass, terms = c("gravity [all]", "GSA")) %>%
#   as.data.frame() %>% 
#   mutate(metric = "Gravity")

# =========================
# 5. COMBINE and clean column names
# =========================
pred_all <- bind_rows(pred_sst, pred_so, pred_so_sd, pred_chla, pred_sbt, pred_depth, pred_fpi_otb) %>%
  rename(value = x, biomass = predicted, lower = conf.low, upper = conf.high, GSA = group) %>%
  # order metrics
  mutate(metric = factor(metric, levels = c("SST mean", "Salinity mean", "Salinity variability", "SBT mean", "Chlorophyll-a mean", "Depth", "Fishing pressure (OTB)")))

# Extract slopes for each metric
slopes_all <- as.data.frame(marginaleffects::avg_slopes(m.biomass, variables = c("sst_mean", "so_mean", "so_sd", "chla_mean", "sbt_mean", "depth", "FPI_OTB"), by = "GSA"))

metric_display_map <- c(
  "sst_mean" = "SST mean",
  #"sst_sd" = "SST variability",
  "so_mean" = "Salinity mean",
  "so_sd" = "Salinity variability",
  "chla_mean" = "Chlorophyll-a mean",
  "sbt_mean" = "SBT mean",
  "depth" = "Depth",
  "FPI_OTB" = "Fishing pressure (OTB)")
  #"FPI_SSF" = "Fishing pressure (SSF)",
  #"gravity" = "Gravity")

gsa_display_map <- c(
  "Catalan_Sea"  = "GSA 6",
  "Gulf_of_Lion" = "GSA 7",
  "Adriatic_Sea" = "GSAs 17-18",
  "Aegean_Sea"   = "GSAs 22-23"
)

# Order controls bottom-to-top stacking (bottom = line 0)
gsa_stack_order <- c("Aegean_Sea", "Adriatic_Sea", "Gulf_of_Lion", "Catalan_Sea")


label_df <- slopes_all %>%
  group_by(term) %>%
  slice(1) %>%
  dplyr::mutate(term_str = as.character(term)) %>%
  dplyr::mutate(
    metric = metric_display_map[term_str],
    stars = dplyr::case_when(
      p.value <= 0.001 ~ "***",
      p.value <= 0.01  ~ "**",
      p.value <= 0.05  ~ "*",
      p.value <= 0.1 ~ ".",
      TRUE ~ "ns"
    ),
    label = paste0("\u03b2 = ", signif(estimate, 2), " ", stars)) %>%
  dplyr::mutate(metric = factor(metric, levels = c("SST mean", "Salinity mean", "Salinity variability", "SBT mean", "Chlorophyll-a mean", "Depth", "Fishing pressure (OTB)")))

not_retained_label <- data.frame(
  metric = factor(c("SST variability", "Fishing pressure (SSF)", "Gravity"), levels = levels(data_long$metric)),
  x = rep(Inf, 3), y = rep(-2.5, 3), label = rep("Not retained", 3))

# =========================
# 7. MULTI-PANEL FIGURE
# =========================
GSA_colors <- c("Adriatic_Sea" = "#0072B2", "Aegean_Sea" = "#D55E00", "Catalan_Sea" = "#009E73", "Gulf_of_Lion" = "#CC79A7")

figure_S4_2 <- ggplot(pred_all, aes(x = value, y = biomass, color = GSA, fill = GSA)) +
   geom_point(data = data_long, aes(x = value, y = log_biomass, color = GSA),  inherit.aes = FALSE, alpha = 0.08, size = 0.5) +
# geom_hex(data = data_long, aes(x = value, y = log_biomass), inherit.aes = FALSE, bins = 40, alpha = 0.4) +
  geom_line(linewidth = 1.3) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = NA) +
  geom_text(data = label_df, aes(x = Inf, y = -2.5, label = label, color = "black"), hjust = 1.04, size = 3.5, inherit.aes = FALSE, show.legend = FALSE) +
  geom_text(data = not_retained_label, aes(x = x, y = y, label = label, color = "black"), inherit.aes = FALSE, size = 3.5, hjust = 1.04, show.legend = FALSE) +
  facet_wrap(~ metric, scales = "free_x", ncol = 2) +
  labs(x = "Standardized metric value", y = "Predicted standardized log biomass", color = "GSA", fill = "GSA") +
  theme_classic(base_size = 12) +
  theme(strip.text = element_text(face = "bold", size = 13),
    strip.background = element_blank(),
    legend.position = "bottom",
    panel.grid = element_blank(),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    plot.title = element_text(face = "bold", size = 16)) +
    scale_color_manual(values = GSA_colors, breaks = c("Catalan_Sea", "Gulf_of_Lion", "Adriatic_Sea", "Aegean_Sea"), labels = c("Northern Spain\n(GSA 6)", "Gulf of Lion\n(GSA 7)", "Adriatic Sea\n(GSAs 17-18)", "Aegean Sea\n(GSAs 22-23)")) +
    scale_fill_manual(values = GSA_colors, breaks = c("Catalan_Sea", "Gulf_of_Lion", "Adriatic_Sea", "Aegean_Sea"), labels = c("Northern Spain\n(GSA 6)", "Gulf of Lion\n(GSA 7)", "Adriatic Sea\n(GSAs 17-18)", "Aegean Sea\n(GSAs 22-23)"))

figure_S4_2
# =========================
# 8. EXPORT
# =========================
ggsave("Figures/Figure_S5.png", plot = figure_S4_2, width = 7, height = 9, dpi = 300)

