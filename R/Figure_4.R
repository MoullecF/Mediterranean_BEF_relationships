# ==============================================================================
# Figure 4: Biomass-biodiversity relationships
#
#   Refits the biomass model and plots predicted biomass against biodiversity metrics, with slopes by GSA.
#
# ==============================================================================

# ==============================================================================
# 1. Load libraries and prepare the environment
# ==============================================================================

source(here::here("R/0-Load libraries.R"))

rm(list = ls())

# ==============================================================================
# 2. Load and prepare data
# ==============================================================================

load(here::here("Data/SEM_data_all_GSA.RData"))

# =========================
# 2. PREPARE DATA
# =========================
SEM_data <- SEM_data %>%
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
  select(log_biomass, GSA, feve, fric, tric, teve, log_Length_max_cwm) %>%
  tidyr::pivot_longer(
    cols = c(feve, fric, tric, teve, log_Length_max_cwm),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(metric = recode(metric,
    feve = "Functional evenness",
    fric = "Functional richness",
    tric = "Taxonomic richness",
    teve = "Taxonomic evenness",
    log_Length_max_cwm = "Log CWM length"
  )) %>%
  mutate(metric = factor(metric, levels = c("Taxonomic richness", "Taxonomic evenness", "Functional richness", "Functional evenness", "Log CWM length")))

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

pred_feve <- ggpredict(m.biomass, terms = c("feve [all]", "GSA")) %>%
  as.data.frame() %>%
  mutate(metric = "Functional evenness")

pred_fric <- ggpredict(m.biomass, terms = c("fric [all]", "GSA")) %>%
  as.data.frame() %>%
  mutate(metric = "Functional richness")

pred_tric <- ggpredict(m.biomass, terms = c("tric [all]", "GSA")) %>%
  as.data.frame() %>%
  mutate(metric = "Taxonomic richness")

pred_teve <- ggpredict(m.biomass, terms = c("teve [all]", "GSA")) %>%
  as.data.frame() %>%
  mutate(metric = "Taxonomic evenness")

pred_size <- ggpredict(m.biomass, terms = c("log_Length_max_cwm [all]", "GSA")) %>%
  as.data.frame() %>%
  mutate(metric = "Log CWM length")

# =========================
# 5. COMBINE
# =========================
pred_all <- bind_rows(
  pred_feve,
  pred_fric,
  pred_tric,
  pred_teve,
  pred_size)

# =========================
# 6. CLEAN COLUMN NAMES
# =========================
pred_all <- pred_all %>%
  rename(value = x, biomass = predicted, lower = conf.low, upper = conf.high, GSA = group) %>%
  # order metrics
  mutate(metric = factor(metric, levels = c("Taxonomic richness", "Taxonomic evenness", "Functional richness", "Functional evenness", "Log CWM length")))

# Extract slopes for all metrics
slopes_all <- as.data.frame(marginaleffects::avg_slopes(m.biomass, variables = c("feve", "fric", "tric", "teve", "log_Length_max_cwm"), by = "GSA"))

# 2. Map original metrics to display names
metric_display_map <- c(
  "feve" = "Functional evenness",
  "fric" = "Functional richness",
  "tric" = "Taxonomic richness",
  "teve" = "Taxonomic evenness",
  "log_Length_max_cwm" = "Log CWM length")

gsa_display_map <- c(
  "Catalan_Sea"  = "GSA 6",
  "Gulf_of_Lion" = "GSA 7",
  "Adriatic_Sea" = "GSAs 17-18",
  "Aegean_Sea"   = "GSAs 22-23")

gsa_stack_order <- c("Aegean_Sea", "Adriatic_Sea", "Gulf_of_Lion", "Catalan_Sea")

# 3. Clean up labels and set up vertical text stacking
label_df <- slopes_all %>%
  dplyr::mutate(
    # Force 'term' and 'GSA' to strings so string operations succeed
    term_str = as.character(term),
    gsa_str  = as.character(GSA)
  ) %>%
  dplyr::mutate(
    metric = metric_display_map[term_str],
    gsa_display = gsa_display_map[gsa_str],
    stars = dplyr::case_when(
      p.value <= 0.001 ~ "***",
      p.value <= 0.01  ~ "**",
      p.value <= 0.05  ~ "*",
      p.value <= 0.1   ~ ".",
      TRUE ~ "ns"
    ),
    label = paste0("\u03b2 (", gsa_display, ") = ", signif(estimate, 2), " ", stars),
    line_num = match(gsa_str, gsa_stack_order) - 1L,
    vjust_val = -(line_num * 1.35 + 0.5)
  ) %>%
  dplyr::mutate(metric = factor(metric, levels = c("Taxonomic richness", "Taxonomic evenness", "Functional richness", "Functional evenness", "Log CWM length")))


# =========================
# 7. MULTI-PANEL FIGURE
# =========================
GSA_colors <- c("Adriatic_Sea" = "#0072B2", "Aegean_Sea" = "#D55E00", "Catalan_Sea" = "#009E73", "Gulf_of_Lion" = "#CC79A7")

figure_4 <- ggplot(pred_all, aes(x = value, y = biomass, color = GSA, fill = GSA)) +
   geom_point(data = data_long, aes(x = value, y = log_biomass, color = GSA),  inherit.aes = FALSE, alpha = 0.08, size = 0.5) +
  geom_line(linewidth = 1.3) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = NA) +
  geom_text(data = label_df, aes(x = Inf, y = -Inf, label = label, color = GSA, vjust = vjust_val), hjust = 1.04, size = 2.5, inherit.aes = FALSE, show.legend = FALSE) +
  facet_wrap(~ metric, scales = "free_x", ncol = 2) +
  labs(x = "Standardized biodiversity metric value", y = "Predicted standardized log biomass", color = "GSA", fill = "GSA") +
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

# =========================
# 8. EXPORT
# =========================
# ggsave(here::here("Figures/Figure_4.png"), plot = figure_4, width = 8, height = 8, dpi = 300)

