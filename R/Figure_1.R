# ==============================================================================
# Script: Figure 1: Trawl positions by GSA
#
#   Maps the positions of the trawl hauls in each GSA
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

# Define GSA colors
GSA_colors <- c("Adriatic_Sea" = "#0072B2", "Aegean_Sea" = "#D55E00", "Catalan_Sea" = "#009E73", "Gulf_of_Lion" = "#CC79A7")

# Count number of hauls per GSA
SEM_data %>%
  group_by(GSA) %>%
  summarise(n = n())

# Create a data frame with the total number of hauls in each GSA to add as text on the map
haul_labels <- data.frame(
  longitude = c(2.3, 5.1, 18.5, 25.2),
  latitude  = c(40.2, 42.3, 39.5, 36),
  label = c(
    "n = 1888",
    "n = 1417",
    "n = 4748",
    "n = 1981"), GSA = c("Catalan_Sea", "Gulf_of_Lion", "Adriatic_Sea", "Aegean_Sea"))

# Plot trawling positions in each GSA in the Mediterranean Sea
 gg_trawling_positions <- ggplot(SEM_data, aes(x = longitude, y = latitude, color = GSA)) +
  geom_jitter(data = SEM_data, aes(x = longitude, y = latitude), width = 0.1, height = 0.1, size = 0.1) +
  scale_color_manual(values = GSA_colors, breaks = c("Catalan_Sea", "Gulf_of_Lion", "Adriatic_Sea", "Aegean_Sea"), labels = c("Northern Spain\n(GSA 6)", "Gulf of Lion\n(GSA 7)", "Adriatic Sea\n(GSAs 17-18)", "Aegean Sea\n(GSAs 22-23)")) +
  guides(color = guide_legend(override.aes = list(size = 3.5))) +
  # Add total number of hauls in each GSA in text on the map
  geom_text(data = haul_labels, aes(x = longitude, y = latitude, label = label, color = GSA), size = 2.5, fontface = "bold", inherit.aes = FALSE) +
  geom_sf(data = rnaturalearth::ne_countries(scale = "large", returnclass = "sf"), fill = "grey90", color = "grey20", inherit.aes = FALSE) +
  scale_x_continuous(limits = range(SEM_data$longitude, na.rm = TRUE) + c(-0.5, 0.5)) +
  scale_y_continuous(limits = range(SEM_data$latitude, na.rm = TRUE) + c(-0.5, 0.5)) +
  coord_sf(expand = FALSE) +
  theme_bw() +
  labs(x = "Longitude", y = "Latitude") + # title = "Trawling positions in each GSA in the Mediterranean Sea", 
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(size = 12),
        plot.title = element_text(size = 12, face = "bold"),
        axis.title = element_text(size = 10))

# ggsave(here::here("Figures/Figure_1.png"), plot = gg_trawling_positions, width = 8, height = 6, dpi = 300)
