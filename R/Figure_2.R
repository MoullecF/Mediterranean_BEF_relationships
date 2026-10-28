# ==============================================================================
# Figure 2: Maps of biodiversity and biomass metrics
#
#   Maps the biodiversity and biomass metrics at haul level
#
# ==============================================================================

# ==============================================================================
# 1. Load libraries and prepare the environment
# ==============================================================================

source("./0-Load libraries.R")

rm(list = ls())

# ==============================================================================
# 2. Load and prepare data
# ==============================================================================

load("Data/SEM_data_all_GSA.RData")

SEM_data$GSA <- factor(SEM_data$GSA, levels = c("Adriatic_Sea", "Aegean_Sea", "Catalan_Sea", "Gulf_of_Lion"))

SEM_data <- SEM_data %>%
 # recode variable GSA to have more readable names
  dplyr::mutate(GSA = recode(GSA,
                              "Adriatic_Sea" = "Adriatic Sea",
                              "Aegean_Sea"   = "Aegean Sea",
                              "Catalan_Sea"  = "Northern Spain",
                              "Gulf_of_Lion" = "Gulf of Lion")) %>%
  # rename columns to have more readable names
  dplyr::rename(
    Latitude = latitude,
    Longitude = longitude,
    Log_biomass = log_biomass,
    Taxonomic_richness = tric,
    Taxonomic_evenness = teve,
    Functional_richness = fric,
    Functional_evenness = feve,
    CWM_length = Length_max_cwm
  ) %>%
  # log transformed CWM_length
  dplyr::mutate(Log_CWM_length = log(CWM_length))

plot_biomass <- ggplot()+
  geom_jitter(data = SEM_data, aes(x = Longitude, y = Latitude, color = Log_biomass), width = 0.1, height = 0.1, size = 0.4, shape = 21) +
  geom_sf(data = rnaturalearth::ne_countries(scale = "large", returnclass = "sf"), fill = "grey90", color = "grey20", inherit.aes = FALSE) +
  scale_color_scico(palette = "vik",
                  midpoint = median(SEM_data$Log_biomass, na.rm = TRUE),
                  name = "Log Biomass") +
  scale_x_continuous(limits = range(SEM_data$Longitude, na.rm = TRUE)+c(-0.5, 0.5)) +
  scale_y_continuous(limits = range(SEM_data$Latitude, na.rm = TRUE)+c(-0.5, 0.5)) +
  coord_sf(expand = FALSE) +
  theme_bw(base_size = 12) +
  theme(legend.position = c(1, 1),   
    legend.direction = "horizontal",
    legend.justification = c("right", "top"),
    legend.background = element_rect(fill = "white", color = "grey50", linewidth = 0.3),
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 12),
    panel.spacing = unit(0, "cm"),
   axis.title.x = element_blank())+
  guides(color = guide_colorbar(title.position = "top", title.hjust = 0.5, barwidth = unit(5, "cm"), barheight = unit(0.4, "cm"))) 

plot_tric <- ggplot()+
  geom_jitter(data = SEM_data, aes(x = Longitude, y = Latitude, color = Taxonomic_richness), width = 0.1, height = 0.1, size = 0.4, shape = 21) +
  geom_sf(data = rnaturalearth::ne_countries(scale = "large", returnclass = "sf"), fill = "grey90", color = "grey20", inherit.aes = FALSE) +
  scale_color_scico(palette = "roma",   
                  midpoint = median(SEM_data$Taxonomic_richness, na.rm = TRUE),
                  name = "Taxonomic richness", direction = -1) +
  scale_x_continuous(limits = range(SEM_data$Longitude, na.rm = TRUE)+c(-0.5, 0.5)) +
  scale_y_continuous(limits = range(SEM_data$Latitude, na.rm = TRUE)+c(-0.5, 0.5)) +
  coord_sf(expand = FALSE) +
  theme_bw(base_size = 12) +
  theme(legend.position = c(1, 1),   
    legend.direction = "horizontal",
    legend.justification = c("right", "top"),
    legend.background = element_rect(fill = "white", color = "grey50", linewidth = 0.3),
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 12),
    panel.spacing = unit(0, "cm"),
    axis.title.y = element_blank(),
    axis.title.x = element_blank())+
  guides(color = guide_colorbar(title.position = "top", title.hjust = 0.5, barwidth = unit(5, "cm"), barheight = unit(0.4, "cm"))) 

plot_teve <- ggplot()+
  geom_jitter(data = SEM_data, aes(x = Longitude, y = Latitude, color = Taxonomic_evenness), width = 0.1, height = 0.1, size = 0.4, shape = 21) +
  geom_sf(data = rnaturalearth::ne_countries(scale = "large", returnclass = "sf"), fill = "grey90", color = "grey20", inherit.aes = FALSE) +
  scale_color_scico(palette = "roma",   
                  midpoint = median(SEM_data$Taxonomic_evenness, na.rm = TRUE),
                  name = "Taxonomic evenness", direction = -1) +
  scale_x_continuous(limits = range(SEM_data$Longitude, na.rm = TRUE)+c(-0.5, 0.5)) +
  scale_y_continuous(limits = range(SEM_data$Latitude, na.rm = TRUE)+c(-0.5, 0.5)) +
  coord_sf(expand = FALSE) +
  theme_bw(base_size = 12) +
  theme(legend.position = c(1, 1),   
    legend.direction = "horizontal",
    legend.justification = c("right", "top"),
    legend.background = element_rect(fill = "white", color = "grey50", linewidth = 0.3),
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 12),
    panel.spacing = unit(0, "cm"),
    axis.title.x = element_blank())+
  guides(color = guide_colorbar(title.position = "top", title.hjust = 0.5, barwidth = unit(5, "cm"), barheight = unit(0.4, "cm")))

plot_fric <- ggplot()+
  geom_jitter(data = SEM_data, aes(x = Longitude, y = Latitude, color = Functional_richness), width = 0.1, height = 0.1, size = 0.4, shape = 21) +
  geom_sf(data = rnaturalearth::ne_countries(scale = "large", returnclass = "sf"), fill = "grey90", color = "grey20", inherit.aes = FALSE) +
  scale_color_scico(palette = "roma",   
                  midpoint = median(SEM_data$Functional_richness, na.rm = TRUE),
                  name = "Functional richness", direction = -1) +
  scale_x_continuous(limits = range(SEM_data$Longitude, na.rm = TRUE)+c(-0.5, 0.5)) +
  scale_y_continuous(limits = range(SEM_data$Latitude, na.rm = TRUE)+c(-0.5, 0.5)) +
  coord_sf(expand = FALSE) +
  theme_bw(base_size = 12) +
  theme(legend.position = c(1, 1),   
    legend.direction = "horizontal",
    legend.justification = c("right", "top"),
    legend.background = element_rect(fill = "white", color = "grey50", linewidth = 0.3),
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 12),
    panel.spacing = unit(0, "cm"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank())+
  guides(color = guide_colorbar(title.position = "top", title.hjust = 0.5, barwidth = unit(5, "cm"), barheight = unit(0.4, "cm")))

plot_feve <- ggplot()+
  geom_jitter(data = SEM_data, aes(x = Longitude, y = Latitude, color = Functional_evenness), width = 0.1, height = 0.1, size = 0.4, shape = 21) +
  geom_sf(data = rnaturalearth::ne_countries(scale = "large", returnclass = "sf"), fill = "grey90", color = "grey20", inherit.aes = FALSE) +
  scale_color_scico(palette = "roma",   
                  midpoint = median(SEM_data$Functional_evenness, na.rm = TRUE),
                  name = "Functional evenness", direction = -1) +
  scale_x_continuous(limits = range(SEM_data$Longitude, na.rm = TRUE)+c(-0.5, 0.5)) +
  scale_y_continuous(limits = range(SEM_data$Latitude, na.rm = TRUE)+c(-0.5, 0.5)) +
  coord_sf(expand = FALSE) +
  theme_bw(base_size = 12) +
  theme(legend.position = c(1, 1),   
    legend.direction = "horizontal",
    legend.justification = c("right", "top"),
    legend.background = element_rect(fill = "white", color = "grey50", linewidth = 0.3),
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 12),
    panel.spacing = unit(0, "cm"))+
  guides(color = guide_colorbar(title.position = "top", title.hjust = 0.5, barwidth = unit(5, "cm"), barheight = unit(0.4, "cm")))

plot_CWM_length <- ggplot()+
  geom_jitter(data = SEM_data, aes(x = Longitude, y = Latitude, color = Log_CWM_length), width = 0.1, height = 0.1, size = 0.4, shape = 21) +
  geom_sf(data = rnaturalearth::ne_countries(scale = "large", returnclass = "sf"), fill = "grey90", color = "grey20", inherit.aes = FALSE) +
  scale_color_scico(palette = "roma",   
                  midpoint = median(SEM_data$Log_CWM_length, na.rm = TRUE),
                  name = "Log CWM length", direction = -1) +
  scale_x_continuous(limits = range(SEM_data$Longitude, na.rm = TRUE)+c(-0.5, 0.5)) +
  scale_y_continuous(limits = range(SEM_data$Latitude, na.rm = TRUE)+c(-0.5, 0.5)) +
  coord_sf(expand = FALSE) +
  theme_bw(base_size = 12) +
  theme(legend.position = c(1, 1),   
    legend.direction = "horizontal",
    legend.justification = c("right", "top"),
    legend.background = element_rect(fill = "white", color = "grey50", linewidth = 0.3),
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 12),
    panel.spacing = unit(0, "cm"),
    axis.title.y = element_blank())+
  guides(color = guide_colorbar(title.position = "top", title.hjust = 0.5, barwidth = unit(5, "cm"), barheight = unit(0.4, "cm")))

# Merge plots into a single figure
combined_plot <- (plot_biomass + plot_tric) /
                 (plot_teve   + plot_fric)  /
                 (plot_feve   + plot_CWM_length) &
  theme(plot.margin = margin(0, 0, 0, 0)) 

combined_plot <- combined_plot +
  plot_layout(axes = "collect", axis_titles = "collect") +  
  plot_annotation(tag_levels = "A")

# ggsave("Figures/Figure_2.png", plot = combined_plot, width = 16, height = 11, dpi = 300)
