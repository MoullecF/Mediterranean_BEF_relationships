library(rstatix)
library(multcompView)
library(ggplot2)
library(dplyr)

load("Data/SEM_data_all_GSA.RData")

# drop NAs if needed
options(na.action = "na.omit")

SEM_data <- SEM_data %>%
  dplyr::mutate(log_Length_max_cwm = log(Length_max_cwm)) %>%
  dplyr::mutate(GSA = dplyr::recode(GSA, "Adriatic_Sea" = "Adriatic Sea", "Aegean_Sea" = "Aegean Sea", "Catalan_Sea" = "Northern Spain", "Gulf_of_Lion" = "Gulf of Lion")) %>%
  dplyr::mutate(GSA = factor(GSA, levels = c("Northern Spain", "Gulf of Lion", "Adriatic Sea", "Aegean Sea"))) %>%
  dplyr::rename(
    Latitude = latitude,
    Longitude = longitude,
    Log_biomass = log_biomass,
    Taxonomic_richness = tric,
    Taxonomic_evenness = teve,
    Functional_richness = fric,
    Functional_evenness = feve,
    CWM_length = Length_max_cwm,
    Log_CWM_length = log_Length_max_cwm)

# Compute mean, standard deviation and median of each metric by GSA
metrics_summary <- SEM_data %>%
  group_by(GSA) %>%
  summarise(
    mean_Log_biomass = round(mean(Log_biomass, na.rm = TRUE), 2),
    sd_Log_biomass = round(sd(Log_biomass, na.rm = TRUE), 2),
    median_Log_biomass = round(median(Log_biomass, na.rm = TRUE), 2),
    mean_Taxonomic_richness = round(mean(Taxonomic_richness, na.rm = TRUE), 2),
    sd_Taxonomic_richness = round(sd(Taxonomic_richness, na.rm = TRUE), 2),
    median_Taxonomic_richness = round(median(Taxonomic_richness, na.rm = TRUE), 2),
    mean_Taxonomic_evenness = round(mean(Taxonomic_evenness, na.rm = TRUE), 2),
    sd_Taxonomic_evenness = round(sd(Taxonomic_evenness, na.rm = TRUE), 2),
    median_Taxonomic_evenness = round(median(Taxonomic_evenness, na.rm = TRUE), 2),
    mean_Functional_richness = round(mean(Functional_richness, na.rm = TRUE), 2),
    sd_Functional_richness = round(sd(Functional_richness, na.rm = TRUE), 2),
    median_Functional_richness = round(median(Functional_richness, na.rm = TRUE), 2),
    mean_Functional_evenness = round(mean(Functional_evenness, na.rm = TRUE), 2),
    sd_Functional_evenness = round(sd(Functional_evenness, na.rm = TRUE), 2),
    median_Functional_evenness = round(median(Functional_evenness, na.rm = TRUE), 2),
    mean_CWM_length = round(mean(CWM_length, na.rm = TRUE), 2),
    sd_CWM_length = round(sd(CWM_length, na.rm = TRUE), 2),
    median_CWM_length = round(median(CWM_length, na.rm = TRUE), 2),
    mean_Log_CWM_length = round(mean(Log_CWM_length, na.rm = TRUE), 2),
    sd_Log_CWM_length = round(sd(Log_CWM_length, na.rm = TRUE), 2),
    median_Log_CWM_length = round(median(Log_CWM_length, na.rm = TRUE), 2)
  )

# Compute statistical differences between GSAs for each metric
anova_results <- list()
metrics <- c("Log_biomass", "Taxonomic_richness", "Taxonomic_evenness", "Functional_richness", "Functional_evenness", "Log_CWM_length")

for (metric in metrics) {
  formula <- as.formula(paste(metric, "~ GSA"))
  anova_results[[metric]] <- oneway.test(formula, data = SEM_data, var.equal = FALSE)
}

# pairwise comparisons for each metric
# Because the assumption of homogeneity of variances was violated (Levene's test, p < 0.001),
# differences among GSA groups were tested using Welch's ANOVA followed by Games–Howell post hoc pairwise comparisons
pairwise_results <- list()
for (metric in metrics) {
  formula <- as.formula(paste(metric, "~ GSA"))
  pairwise_results[[metric]] <- games_howell_test(SEM_data, formula)
}

## Function to get letters for pairwise comparisons
# Build named p-value vector from pairwise results
pw_logBiomass <- pairwise_results$Log_biomass
pw_taxoR <- pairwise_results$Taxonomic_richness
pw_taxoEv <- pairwise_results$Taxonomic_evenness
pw_funcR <- pairwise_results$Functional_richness
pw_funcEv <- pairwise_results$Functional_evenness
pw_CWM <- pairwise_results$Log_CWM_length

p_vec_logBiomass <- setNames(pw_logBiomass$p.adj, paste(pw_logBiomass$group1, pw_logBiomass$group2, sep = "-"))
p_vec_taxoR <- setNames(pw_taxoR$p.adj, paste(pw_taxoR$group1, pw_taxoR$group2, sep = "-"))
p_vec_taxoEv <- setNames(pw_taxoEv$p.adj, paste(pw_taxoEv$group1, pw_taxoEv$group2, sep = "-"))
p_vec_funcR <- setNames(pw_funcR$p.adj, paste(pw_funcR$group1, pw_funcR$group2, sep = "-"))
p_vec_funcEv <- setNames(pw_funcEv$p.adj, paste(pw_funcEv$group1, pw_funcEv$group2, sep = "-"))
p_vec_CWM <- setNames(pw_CWM$p.adj, paste(pw_CWM$group1, pw_CWM$group2, sep = "-"))

# Compute compact letter display (CLD)
cld_logBiomass <- multcompLetters(p_vec_logBiomass)$Letters
cld_taxoR <- multcompLetters(p_vec_taxoR)$Letters
cld_taxoEv <- multcompLetters(p_vec_taxoEv)$Letters
cld_funcR <- multcompLetters(p_vec_funcR)$Letters
cld_funcEv <- multcompLetters(p_vec_funcEv)$Letters
cld_CWM <- multcompLetters(p_vec_CWM)$Letters

# Build annotation df with y-position placed above each error bar
letters_df_logBiomass <- tibble(GSA = names(cld_logBiomass), letter = cld_logBiomass) %>%
  left_join(metrics_summary %>%
      mutate(GSA = as.character(GSA)) %>% 
      select(GSA, mean_Log_biomass, sd_Log_biomass),
    by = "GSA") %>%
  mutate(y_pos = mean_Log_biomass + sd_Log_biomass + 0.05*(mean_Log_biomass + sd_Log_biomass))

letters_df_taxoR <- tibble(GSA = names(cld_taxoR), letter = cld_taxoR) %>%
  left_join(metrics_summary %>%
      mutate(GSA = as.character(GSA)) %>% 
      select(GSA, mean_Taxonomic_richness, sd_Taxonomic_richness),
    by = "GSA") %>%
  mutate(y_pos = mean_Taxonomic_richness + sd_Taxonomic_richness + 0.05*(mean_Taxonomic_richness + sd_Taxonomic_richness))

letters_df_taxoEv <- tibble(GSA = names(cld_taxoEv), letter = cld_taxoEv) %>%
  left_join(metrics_summary %>%
      mutate(GSA = as.character(GSA)) %>% 
      select(GSA, mean_Taxonomic_evenness, sd_Taxonomic_evenness),
    by = "GSA") %>%
  mutate(y_pos = mean_Taxonomic_evenness + sd_Taxonomic_evenness + 0.05*(mean_Taxonomic_evenness + sd_Taxonomic_evenness))

letters_df_funcR <- tibble(GSA = names(cld_funcR), letter = cld_funcR) %>%
  left_join(metrics_summary %>%
      mutate(GSA = as.character(GSA)) %>% 
      select(GSA, mean_Functional_richness, sd_Functional_richness),
    by = "GSA") %>%
  mutate(y_pos = mean_Functional_richness + sd_Functional_richness + 0.05*(mean_Functional_richness + sd_Functional_richness))

letters_df_funcEv <- tibble(GSA = names(cld_funcEv), letter = cld_funcEv) %>%
  left_join(metrics_summary %>%
      mutate(GSA = as.character(GSA)) %>% 
      select(GSA, mean_Functional_evenness, sd_Functional_evenness),
    by = "GSA") %>%
  mutate(y_pos = mean_Functional_evenness + sd_Functional_evenness + 0.05*(mean_Functional_evenness + sd_Functional_evenness))

letters_df_CWM <- tibble(GSA = names(cld_CWM), letter = cld_CWM) %>%
  left_join(metrics_summary %>%
      mutate(GSA = as.character(GSA)) %>% 
      select(GSA, mean_Log_CWM_length, sd_Log_CWM_length),
    by = "GSA") %>%
  mutate(y_pos = mean_Log_CWM_length + sd_Log_CWM_length + 0.05*(mean_Log_CWM_length + sd_Log_CWM_length))

# Plot
GSA_colors = c("Adriatic Sea" = "#0072B2","Aegean Sea" = "#D55E00", "Northern Spain" = "#009E73", "Gulf of Lion" = "#CC79A7")

Log_biomass_gsa <- ggplot(metrics_summary, aes(x = GSA, y = mean_Log_biomass, fill = GSA)) +
  geom_col(width = 0.7) +
  geom_errorbar(aes(ymin = mean_Log_biomass - sd_Log_biomass, ymax = mean_Log_biomass + sd_Log_biomass), width = 0.1) +
  geom_text(data = letters_df_logBiomass,
    aes(x = GSA, y = y_pos, label = letter),
    inherit.aes = FALSE,
    size = 4.5,
    fontface = "bold") +
  scale_fill_manual(values = GSA_colors) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic() +
  labs(x = NULL, y = "Log biomass (g/ha)") +
  theme(legend.position = "none")
Log_biomass_gsa

taxo_rich_gsa <- ggplot(metrics_summary, aes(x = GSA, y = mean_Taxonomic_richness, fill = GSA)) +
  geom_col(width = 0.7) +
  geom_errorbar(aes(ymin = mean_Taxonomic_richness - sd_Taxonomic_richness, ymax = mean_Taxonomic_richness + sd_Taxonomic_richness), width = 0.1) +
  geom_text(data = letters_df_taxoR,
    aes(x = GSA, y = y_pos, label = letter),
    inherit.aes = FALSE,
    size = 4.5,
    fontface = "bold") +
  scale_fill_manual(values = GSA_colors) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic() +
  labs(x = NULL, y = "Taxonomic richness") +
  theme(legend.position = "none")
taxo_rich_gsa

taxo_ev_gsa <- ggplot(metrics_summary, aes(x = GSA, y = mean_Taxonomic_evenness, fill = GSA)) +
  geom_col(width = 0.7) +
  geom_errorbar(aes(ymin = mean_Taxonomic_evenness - sd_Taxonomic_evenness, ymax = mean_Taxonomic_evenness + sd_Taxonomic_evenness), width = 0.1) +
  geom_text(data = letters_df_taxoEv,
    aes(x = GSA, y = y_pos, label = letter),
    inherit.aes = FALSE,
    size = 4.5,
    fontface = "bold") +
  scale_fill_manual(values = GSA_colors) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic(base_size = 12) +
  labs(x = NULL, y = "Taxonomic evenness") +
  theme(legend.position = "none")
taxo_ev_gsa

Func_rich_gsa <- ggplot(metrics_summary, aes(x = GSA, y = mean_Functional_richness, fill = GSA)) +
  geom_col(width = 0.7) +
  geom_errorbar(aes(ymin = mean_Functional_richness - sd_Functional_richness, ymax = mean_Functional_richness + sd_Functional_richness), width = 0.1) +
  geom_text(data = letters_df_funcR,
    aes(x = GSA, y = y_pos, label = letter),
    inherit.aes = FALSE,
    size = 4.5,
    fontface = "bold") +
  scale_fill_manual(values = GSA_colors) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic(base_size = 12) +
  labs(x = NULL, y = "Functional richness") +
  theme(legend.position = "none")
Func_rich_gsa

Func_ev_gsa <- ggplot(metrics_summary, aes(x = GSA, y = mean_Functional_evenness, fill = GSA)) +
  geom_col(width = 0.7) +
  geom_errorbar(aes(ymin = mean_Functional_evenness - sd_Functional_evenness, ymax = mean_Functional_evenness + sd_Functional_evenness), width = 0.1) +
  geom_text(data = letters_df_funcEv,
    aes(x = GSA, y = y_pos, label = letter),
    inherit.aes = FALSE,
    size = 4.5,
    fontface = "bold") +
  scale_fill_manual(values = GSA_colors) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic(base_size = 12) +
  labs(x = NULL, y = "Functional evenness") +
  theme(legend.position = "none")
Func_ev_gsa

CWM_length_gsa <- ggplot(metrics_summary, aes(x = GSA, y = mean_Log_CWM_length, fill = GSA)) +
  geom_col(width = 0.7) +
  geom_errorbar(aes(ymin = mean_Log_CWM_length - sd_Log_CWM_length, ymax = mean_Log_CWM_length + sd_Log_CWM_length), width = 0.1) +
  geom_text(data = letters_df_CWM,
    aes(x = GSA, y = y_pos, label = letter),
    inherit.aes = FALSE,
    size = 4.5,
    fontface = "bold") +
  scale_fill_manual(values = GSA_colors) +
  scale_y_continuous(expand = c(0, 0.1)) +
  theme_classic(base_size = 12) +
  labs(x = NULL, y = "Log CWM length (cm)") +
  theme(legend.position = "none")
CWM_length_gsa

# Combine all plots into a single figure
combined_plot <- (Log_biomass_gsa | CWM_length_gsa) / (taxo_rich_gsa | taxo_ev_gsa) / (Func_rich_gsa | Func_ev_gsa) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A")
combined_plot

ggsave("Figures/Figure_S2.png", plot = combined_plot, width = 12, height = 17, dpi = 300)
