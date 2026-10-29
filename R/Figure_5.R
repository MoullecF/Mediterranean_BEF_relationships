# ==============================================================================
# Figure 5: Direct, indirect and total effects
#
#   Combined figure of direct, indirect and total effects of predictors and biodiversity on biomass.
#
# ==============================================================================

# ==============================================================================
# 1. Load libraries and prepare the environment
# ==============================================================================

source(here::here("R/0-Load libraries.R"))

rm(list = ls())

# Load necessary data and functions

load(here::here("Tools/extract_direct_effects_function.RData"))
load(here::here("Tools/compute_indirect_effects_GSA_decomp_function.RData"))
load(here::here("Tools/compute_total_effects_GSA_function.RData"))

load(here::here("Outputs/SEM_results_ALL_GSA.RData"))
load(here::here("Outputs/SEM_model_ALL_GSA.RData"))
load(here::here("Data/SEM_data_all_GSA.RData"))

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

SEM_data_complete <- SEM_data %>%
  # 1. Create the dummy 0/1 columns
  fastDummies::dummy_cols("GSA", remove_first_dummy = TRUE, remove_selected_columns = FALSE)

# Compute direct effects of biodiversity on log biomass across GSAs
direct_effects_bio <- extract_direct_effects(
  sem_model = SEM_model,
  response = "log_biomass",
  predictors = c("tric", "teve", "feve", "fric", "log_Length_max_cwm"),
  raw_data_object = SEM_data_complete)

# Compute direct effects of environmental and anthropogenic predictors on log biomass across GSAs
direct_effects_env <- extract_direct_effects(
  sem_model = SEM_model,
  response = "log_biomass",
  predictors = c("sst_mean", "sst_sd", "sbt_mean", "so_mean", "so_sd", "chla_mean", "depth", "FPI_OTB", "FPI_SSF", "gravity"),
  raw_data_object = SEM_data_complete)

# Compute indirect effects of predictors on log biomass across GSAs
map_single_mediator <- function(mediator) {
  dplyr::case_when(
    mediator %in% c("fric", "feve")         ~ "Functional diversity",
    mediator %in% c("tric", "teve")         ~ "Taxonomic diversity",
    mediator %in% c("log_Length_max_cwm")       ~ "Traits",
    TRUE                                    ~ "Other")
}

group_fun <- function(mediators_vector) {
  mapped_elements <- map_single_mediator(mediators_vector)
  combined_path_string <- paste(unique(mapped_elements), collapse = " → ")
  dplyr::case_when(
    combined_path_string == "Taxonomic diversity → Functional diversity" ~ "Taxonomic diversity → Functional diversity",
    combined_path_string == "Taxonomic diversity"                       ~ "Taxonomic diversity",
    combined_path_string == "Functional diversity"                      ~ "Functional diversity",
    combined_path_string == "Traits"                                    ~ "Traits",
    TRUE                                                                ~ "Other")
}

indirect_effects <- compute_indirect_effects_GSA_decomp(
  model = SEM_model,
  response = "log_biomass",
  gsa_levels = c("Adriatic_Sea", "Aegean_Sea", "Catalan_Sea", "Gulf_of_Lion"),
  group_fun = group_fun)

# Compute total effects by combining direct and indirect effects
total_effects_bio <- compute_total_effects_GSA(
  direct_df = direct_effects_bio$total_effects,
  indirect_df = indirect_effects
)

total_effects_env <- compute_total_effects_GSA(
  direct_df = direct_effects_env$total_effects,
  indirect_df = indirect_effects
)

# Figure settings
base_text_size  <- 16  # base font size (pt) of all panels: axis titles, axis labels, strips, legend
label_text_size <- 5   # size (mm) of significance labels on bars
tag_text_size   <- 22  # size (pt) of panel letters A-E
GSA_colors = c("Adriatic Sea" = "#0072B2","Aegean Sea" = "#D55E00", "Northern Spain" = "#009E73", "Gulf of Lion" = "#CC79A7")
mediator_colors <- c("Functional diversity" = "#8DA0CB", "Taxonomic diversity" = "#FC8D62", "Traits" = "#66C2A5", "Taxonomic diversity → Functional diversity" = "#E78AC3")

# ==============================================================================
# 1. CLEAN AND FORMAT DATASETS FOR ALL FIVE PANELS
# ==============================================================================

# Unified lookup array for display titles across panels
predictor_clean_map <- c(
  "tric"           = "Taxonomic richness",
  "teve"           = "Taxonomic evenness",
  "fric"           = "Functional richness",
  "feve"           = "Functional evenness",
  "log_Length_max_cwm" = "Log CWM length",
  "sst_mean"       = "SST mean",
  "sst_sd"         = "SST sd",
  "sbt_mean"       = "SBT mean",
  "so_mean"        = "SO mean",
  "so_sd"          = "SO sd",
  "chla_mean"      = "Chla mean",
  "depth"          = "Depth",
  "FPI_OTB"        = "Fishing press. OTB",
  "FPI_SSF"        = "Fishing press. SSF",
  "gravity"        = "Gravity"
)

bio_levels <- c("Taxonomic richness", "Taxonomic evenness", "Functional richness", "Functional evenness", "Log CWM length")
env_levels <- c("SST mean", "SST sd", "SBT mean", "SO mean", "SO sd", "Chla mean", "Depth",
                "Fishing press. OTB", "Fishing press. SSF", "Gravity")

clean_GSA <- function(GSA_level) {
  GSA_clean <- GSA_level %>% stringr::str_remove("^GSA_") %>% stringr::str_replace_all("_", " ")
  GSA_clean <- dplyr::if_else(GSA_clean == "Catalan Sea", "Northern Spain", GSA_clean)
  factor(GSA_clean, levels = c("Northern Spain", "Gulf of Lion", "Adriatic Sea", "Aegean Sea"))
}

# A. DIRECT EFFECTS OF ENVIRONMENTAL AND ANTHROPOGENIC PREDICTORS (ordered by mean effect)
direct_data_env <- direct_effects_env$total_effects %>%
  dplyr::mutate(
    GSA_clean = clean_GSA(GSA_level),
    predictor_clean = predictor_clean_map[predictor]) %>%
  dplyr::group_by(predictor_clean) %>%
  dplyr::mutate(mean_total = mean(total_std, na.rm = TRUE)) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(mean_total) %>%
  dplyr::mutate(predictor_clean = factor(predictor_clean, levels = unique(predictor_clean))) %>%
  dplyr::select(-mean_total) %>%
  # Replace "ns" by "not retained" for gravity, FPI SSF and SST sd not retained by dredge()
  dplyr::mutate(sig_label = dplyr::if_else(predictor_clean %in% c("Gravity", "Fishing press. SSF", "SST sd") & significance_stars == "ns", "not retained", significance_stars))

# B. DIRECT EFFECTS OF BIODIVERSITY
direct_data_bio <- direct_effects_bio$total_effects %>%
  dplyr::mutate(
    GSA_clean = clean_GSA(GSA_level),
    sig_label = dplyr::if_else(!significant, "ns", NA_character_),
    predictor_clean = factor(predictor_clean_map[predictor], levels = bio_levels)
  )

# C. INDIRECT EFFECTS DATA CLEANING (Aggregated by mediator_group for stacked bars)
indirect_data <- indirect_effects$effects %>%
  # Force to character string first so the text pattern matching is 100% reliable
  dplyr::mutate(predictor_str = as.character(predictor)) %>%
  # Cleanly drop the regional baseline dummy variables
  dplyr::filter(!stringr::str_detect(predictor_str, "^GSA_")) %>%
  dplyr::mutate(
    GSA_clean = clean_GSA(GSA_level),
    predictor_clean = factor(predictor_clean_map[predictor_str], levels = rev(c(bio_levels, env_levels)))
  ) %>%
  # Drop the temporary string helper column
  dplyr::select(-predictor_str)

# D. TOTAL EFFECTS OF ENVIRONMENTAL AND ANTHROPOGENIC PREDICTORS
total_data_env <- total_effects_env %>%
  dplyr::mutate(
    GSA_clean = clean_GSA(GSA_level),
    predictor_clean = factor(predictor_clean_map[predictor], levels = env_levels))

# E. TOTAL EFFECTS OF BIODIVERSITY
total_data_bio <- total_effects_bio %>%
  dplyr::mutate(
    GSA_clean = clean_GSA(GSA_level),
    predictor_clean = factor(predictor_clean_map[predictor], levels = bio_levels))

# ==============================================================================
# 2. DESIGN THE FIVE PANELS
# ==============================================================================

# ---- PANEL A: DIRECT EFFECTS OF ENVIRONMENTAL AND ANTHROPOGENIC PREDICTORS (Figure 6, panel A) ----
direct_data_env <- direct_data_env %>%
  mutate(sig_label = dplyr::na_if(sig_label, "."), significance_stars = dplyr::na_if(significance_stars, "."))

plot_A <- ggplot(direct_data_env %>% dplyr::filter(GSA_clean == "Adriatic Sea"), aes(x = predictor_clean, y = total_std)) +
    geom_col(width = 0.7, fill = "#3B9AB2") +
    geom_text(aes(label = dplyr::case_when(is.na(sig_label) ~ significance_stars, TRUE ~ sig_label), hjust = dplyr::if_else(total_std >= 0, -0.3, 1.3)), vjust = .8 , size = label_text_size) +
    geom_point(data = data.frame(predictor_clean = "Chla mean", total_std = 0.034), aes(x = predictor_clean, y = total_std), colour = "black", size = 0.5, inherit.aes = FALSE) +  
    geom_hline(yintercept = 0, linewidth = 0.4, linetype = "dashed", colour = "black") +
    labs(title = "Standardized direct effects of environmental and anthropogenic predictors on log biomass", x = "Predictors", y = "Direct effect") +
    theme_classic(base_size = base_text_size) +
    theme(strip.background  = element_blank(),
      strip.text        = element_text(face = "bold"),
      axis.text.x       = element_text(angle = 0, hjust = 1),
      panel.grid.major.x = element_blank(),
      panel.background   = element_blank(),
      panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3)) +
  coord_flip()

# ---- PANEL B: DIRECT EFFECTS OF BIODIVERSITY (Figure 5, panel A) ----
plot_B <- ggplot(direct_data_bio, aes(x = GSA_clean, y = total_std, fill = GSA_clean)) +
  geom_col(width = 0.7) +
  geom_text(aes(
    label = dplyr::case_when(is.na(sig_label) ~ significance_stars, TRUE ~ sig_label),
    vjust = dplyr::if_else(total_std >= 0, -0.4, 1.4)
  ), size = label_text_size) +
  geom_hline(yintercept = 0, linewidth = 0.4, linetype = "dashed", colour = "black") +
  ylim(-0.45, 0.52) +
  scale_fill_manual(values = GSA_colors, guide = "none") +
  facet_wrap(~ predictor_clean, scales = "free_y", nrow = 1, labeller = label_wrap_gen(width = 25)) +
  labs(title = "Standardized direct effects of biodiversity on log biomass", x = NULL, y = "Direct effect") +
  theme_classic(base_size = base_text_size) +
  theme(
    strip.background   = element_blank(),
    strip.text         = element_text(face = "bold"),
    axis.text.x        = element_blank(), # Remove text to avoid overlap in the stack
    axis.ticks.x       = element_blank(),
    panel.grid.major.y = element_line(color = "grey95", linewidth = 0.3)
  )

# ---- PANEL C: INDIRECT EFFECTS (Figures 5 and 6, panel B) ----
plot_C <- ggplot(indirect_data, aes(x = indirect_std, y = predictor_clean, fill = mediator_group)) +
  geom_col(width = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  scale_fill_manual(values = mediator_colors) +
  facet_wrap(~ GSA_clean, nrow = 1) +
  scale_y_discrete(expand = c(0.05, 0.05)) +
  labs(title = "Standardized indirect effects of predictors on log biomass", x = "Indirect effect", y = "Predictors", fill = "Pathway group") +
  theme_classic(base_size = base_text_size) +
  theme(
    strip.background   = element_blank(),
    strip.text         = element_text(face = "bold"),
    panel.grid.major.y = element_line(color = "grey95", linewidth = 0.3),
    legend.position    = "bottom",
    legend.text        = element_text(size = 18, margin = margin(l= 5,r = 20)),
    legend.title       = element_text(size = 18),
    legend.key.size    = unit(0.8, "cm"),
    legend.spacing.x   = unit(0.8, "cm"))

# ---- PANEL D: TOTAL EFFECTS OF ENVIRONMENTAL AND ANTHROPOGENIC PREDICTORS (Figure 6, panel C) ----
plot_D <- ggplot(total_data_env, aes(x = GSA_clean, y = total_std, fill = GSA_clean)) +
  geom_col(width = 0.7) +
  geom_hline(yintercept = 0, linewidth = 0.4, linetype = "dashed", colour = "black") +
  ylim(-0.33, 0.15) +
  scale_fill_manual(values = GSA_colors, guide = "none") +
  facet_wrap(~ predictor_clean, scales = "free_y", nrow = 1, labeller = label_wrap_gen(width = 12)) +
  labs(title = "Standardized total effects of environmental and anthropogenic predictors on log biomass", x = NULL, y = "Total effect") +
  theme_classic(base_size = base_text_size) +
  theme(
    strip.background   = element_blank(),
    strip.text         = element_text(face = "bold"),
    #axis.text.x        = element_text(angle = 35, hjust = 1),
    axis.text.x        = element_blank(), # Remove text to avoid overlap in the stack
    axis.ticks.x       = element_blank(),
    panel.grid.major.y = element_line(color = "grey95", linewidth = 0.3))

# ---- PANEL E: TOTAL EFFECTS OF BIODIVERSITY (Figure 5, panel C) ----
plot_E <- ggplot(total_data_bio, aes(x = GSA_clean, y = total_std, fill = GSA_clean)) +
  geom_col(width = 0.7) +
  geom_hline(yintercept = 0, linewidth = 0.4, linetype = "dashed", colour = "black") +
  ylim(-0.45, 0.55) +
  scale_fill_manual(values = GSA_colors, guide = "none") +
  facet_wrap(~ predictor_clean, scales = "free_y", nrow = 1, labeller = label_wrap_gen(width = 25)) +
  labs(title = "Standardized total effects of biodiversity on log biomass", x = NULL, y = "Total effect") +
  theme_classic(base_size = base_text_size) +
  theme(
    strip.background   = element_blank(),
    strip.text         = element_text(face = "bold"),
    axis.text.x        = element_text(angle = 35, hjust = 1, size = 18),
    panel.grid.major.y = element_line(color = "grey95", linewidth = 0.3)
  )

# ==============================================================================
# 3. ASSEMBLE COMPLETE MULTI-PANEL FIGURE
# ==============================================================================
# Vertical stack, panels automatically tagged A to E
figure_5_6 <- (plot_A / plot_B / plot_C / plot_D / plot_E) +
  plot_layout(heights = c(1, 0.8, 1.4, 0.8, 0.8)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = tag_text_size))

# Save
# ggsave(here::here("Figures/Figure_5.png"), plot = figure_5_6, width = 18.5, height = 27, dpi = 300)
