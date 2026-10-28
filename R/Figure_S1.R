# ==============================================================================
# Figure S1: Number of hauls per GSA and year
#
#   Plots the number of hauls in each GSA and year.
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

# Plot number of hauls in each GSA and year
SEM_data$id <- 1:nrow(SEM_data)
haul_counts_gsa <- SEM_data %>%
  distinct(year, GSA, id) %>%
  count(year, GSA, name = "haul_count") %>%
  mutate(GSA = factor(GSA, levels = c("Aegean_Sea", "Adriatic_Sea", "Gulf_of_Lion", "Catalan_Sea")))

gg_number_hauls_gsa <- ggplot(haul_counts_gsa, aes(year, GSA, fill = haul_count)) +
  geom_tile(linewidth = 0.1, color = "grey80") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_discrete(labels = c(
      Catalan_Sea = "Northern Spain",
      Gulf_of_Lion = "Gulf of Lion",
      Adriatic_Sea = "Adriatic Sea",
      Aegean_Sea = "Aegean Sea"), expand = c(0, 0)) +
  viridis::scale_fill_viridis(na.value = "white") +
  labs(title = "Number of hauls in each GSA and year", x = "Years", y = NULL, fill = "Number\nof hauls") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.position = "right",
    legend.text = element_text(size = 12),
    plot.title = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text.y = element_text(
      color = GSA_colors[levels(haul_counts_gsa$GSA)],
      face = "bold"))

# ggsave("Figures/Figure_S1.png", plot = gg_number_hauls_gsa, width = 8, height = 4, dpi = 300)
