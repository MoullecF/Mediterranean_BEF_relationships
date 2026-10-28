library(ggplot2)
library(piecewiseSEM)
library(dplyr)
library(tidyr)


rm(list = ls())

# Load SEM model and plot R² values for each sub-model
load("Outputs/SEM_model_ALL_GSA.RData")
  
# Extract R²
r2 <- rsquared(SEM_model)
  
# Prepare data
r2_long <- r2 %>%
  dplyr::mutate(
    Random = Conditional - Marginal, 
    Response = dplyr::case_when( 
      Response == "tric" ~ "Taxonomic\nrichness", 
      Response == "teve" ~ "Taxonomic\nevenness", 
      Response == "fric" ~ "Functional\nrichness", 
      Response == "feve" ~ "Functional\nevenness", 
      Response == "log_Length_max_cwm" ~ "Log\nCWM length", 
      Response == "log_biomass" ~ "Log\nBiomass", 
      TRUE ~ Response)) %>% 
  dplyr::mutate(Response = forcats::fct_reorder(Response, Marginal, .desc = TRUE)) %>% 
  dplyr::select(Response, Marginal, Random) %>% 
  tidyr::pivot_longer(cols = c(Marginal, Random), names_to = "Component", values_to = "R2") %>% 
  dplyr::mutate(Component = factor(Component, levels = c("Marginal", "Random")))
  
# Plot
p <- ggplot(r2_long, aes(x = Response, y = R2, fill = Component)) +
  geom_bar(stat = "identity", color = "black", position = position_stack(reverse = TRUE)) +
  scale_fill_manual(values = c("Marginal" = "#78B7C5", "Random" = "#EBCC2A"),
      labels = c("Marginal" = "Fixed effects", "Random"   = "Random effects")) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 0.67)) +
  labs(x = NULL, y = expression(R^2), fill = NULL) +
  theme_classic(base_size = 13) +
  theme(legend.position = c(0.80, 0.85))

# Save figure
ggsave("Figures/Figure_3_rsquared.png", plot = p, width = 9, height = 4.5, dpi = 300)
