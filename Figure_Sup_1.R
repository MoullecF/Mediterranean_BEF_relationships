# Load libraries
library(ggplot2)
library(dplyr)

rm(list = ls())

# Load data
load("Data/SEM_data_all_GSA.RData")
SEM_data <- SEM_data %>%
  # select only the columns we need for the figure
  select(GSA, log_biomass, tric, teve, fric, feve, Length_max_cwm) %>%
  # pivot longer to have one row per GSA and trait
  tidyr::pivot_longer(cols = c(log_biomass, tric, teve, fric, feve, Length_max_cwm), names_to = "metrics", values_to = "metric_values")

# rename the metrics for better display in the figure
SEM_data$metrics <- recode(SEM_data$metrics,
                           "log_biomass" = "Log Biomass",
                           "tric" = "Taxonomic richness",
                           "teve" = "Taxonomic evenness",
                           "fric" = "Functional richness",
                           "feve" = "Functional evenness",
                           "Length_max_cwm" = "CWM length")

# Rename GSA for better display in the figure
SEM_data$GSA <- recode(SEM_data$GSA,
                       "Adriatic_Sea" = "Adriatic Sea",
                       "Aegean_Sea" = "Aegean Sea",
                       "Catalan_Sea" = "Northern Spain",
                       "Gulf_of_Lion" = "Gulf of Lion")

# reorder the metrics and GSA for better display in the figure
SEM_data$metrics <- factor(SEM_data$metrics, levels = c("Log Biomass", "Taxonomic richness", "Taxonomic evenness", "Functional richness", "Functional evenness", "CWM length"))
SEM_data$GSA <- factor(SEM_data$GSA, levels = c("Northern Spain", "Gulf of Lion", "Adriatic Sea", "Aegean Sea"))


GSA_colors <- c("Adriatic Sea" = "#0072B2", "Aegean Sea" = "#D55E00", "Northern Spain" = "#009E73", "Gulf of Lion" = "#CC79A7")

means <- SEM_data %>%
  group_by(GSA, metrics) %>%
  summarise(mean_metric = mean(metric_values, na.rm = TRUE),
            .groups = "drop")

dens_max <- SEM_data %>%
  group_by(metrics, GSA) %>%
  summarise(mean_metric = mean(metric_values, na.rm = TRUE),
            .groups = "drop") %>%
  left_join(
    SEM_data %>%
      group_by(metrics) %>%
      summarise(
        ymax = {
          d <- density(metric_values, na.rm = TRUE)
          max(d$y)
        },
        .groups = "drop"
      ),
    by = "metrics"
  ) %>%
  mutate(ypos = ymax / 1.5)

gg_density_metric <- ggplot(SEM_data, aes(x = metric_values, fill = GSA)) +
  geom_density(alpha = 0.2, colour = NA) +
  scale_fill_manual(values = GSA_colors) +
  geom_segment(data = dens_max, aes(x = mean_metric, xend = mean_metric, y = 0, yend = ypos), linetype = "dashed", linewidth = 0.5, inherit.aes = FALSE) +
  geom_point(data = dens_max, aes(x = mean_metric, y = ypos, fill = GSA), shape = 21, size = 3, inherit.aes = FALSE)  +
  facet_wrap(~metrics, scales = "free") +
  labs(x = "Metric value", y = "Density")+
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic(base_size = 12) +
  theme(strip.text = element_text(face = "bold", size = 13),
    strip.background = element_blank(),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    panel.grid = element_blank(),
    axis.title = element_text(size = 12),
    axis.text = element_text(color = "black"))

ggsave("Figures/Figure_S5.png", plot = gg_density_metric, width = 14, height = 8, dpi = 300)

# Figure caption: Density plots of the six metrics (log biomass, taxonomic richness, taxonomic evenness, functional richness, functional evenness, and CWM length)
# for the four GSAs (Northern Spain, Gulf of Lion, Adriatic Sea, and Aegean Sea).
# The vertical lines and points indicate the mean value of each metric in each GSA.
