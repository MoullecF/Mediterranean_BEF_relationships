library(rstatix)
library(multcompView)
library(ggplot2)
library(dplyr)

load("Data/SEM_data_all_GSA.RData")

# Other scripts set na.fail globally; formula tests here need to drop NAs (e.g. Biomass_group)
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

# Correlation between biomass and longitude
cor.test(SEM_data$Log_biomass, SEM_data$Longitude, method = "spearman")
ggplot(SEM_data, aes(Longitude, Log_biomass)) +
  geom_point(aes(fill = GSA), shape = 21) +
  geom_smooth(method = "lm", se = FALSE)+
  theme_classic()

# Correlation between biomass and latitude
cor.test(SEM_data$Log_biomass, SEM_data$Latitude, method = "spearman")

# Correlation between biomass and latitude in the Adriatic Sea
SEM_data_Adriatic <- SEM_data %>% filter(GSA == "Adriatic Sea")
cor(SEM_data_Adriatic$Log_biomass, SEM_data_Adriatic$Latitude, use = "complete.obs", method = "spearman")
cor.test(SEM_data_Adriatic$Log_biomass, SEM_data_Adriatic$Latitude, method = "spearman")

# Student's t-test to compare biomass ≥43°N in Adriatic Sea vs. <43°N in Adriatic Sea
SEM_data_Adriatic$Latitude_group <- ifelse(SEM_data_Adriatic$Latitude >= 43, "≥43°N", "<43°N")
round(mean(SEM_data_Adriatic$Log_biomass[SEM_data_Adriatic$Latitude_group == "≥43°N"], na.rm = TRUE), 2)
round(sd(SEM_data_Adriatic$Log_biomass[SEM_data_Adriatic$Latitude_group == "≥43°N"], na.rm = TRUE), 2)
round(mean(SEM_data_Adriatic$Log_biomass[SEM_data_Adriatic$Latitude_group == "<43°N"], na.rm = TRUE), 2)
round(sd(SEM_data_Adriatic$Log_biomass[SEM_data_Adriatic$Latitude_group == "<43°N"], na.rm = TRUE), 2)

car::leveneTest(Log_biomass ~ Latitude_group, data = SEM_data_Adriatic)
t_test_result <- t.test(Log_biomass ~ Latitude_group, data = SEM_data_Adriatic, var.equal = FALSE)
t_test_result

# Correlation between taxonomic richness and Log biomass
cor.test(SEM_data$Taxonomic_richness, SEM_data$Log_biomass, method = "spearman")

# Correlation between taxonomic evenness and Log biomass
cor.test(SEM_data$Taxonomic_evenness, SEM_data$Log_biomass, method = "spearman")
quantile(SEM_data$Log_biomass, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)

# Comparing taxonomic evenness from hauls in the top biomass quartile (≥75th percentile) vs. the bottom biomass quartile (≤25th percentile)
SEM_data$Biomass_group <- ifelse(SEM_data$Log_biomass >= quantile(SEM_data$Log_biomass, 0.75, na.rm = TRUE), "≥75th percentile", 
                                  ifelse(SEM_data$Log_biomass <= quantile(SEM_data$Log_biomass, 0.25, na.rm = TRUE), "≤25th percentile", NA))
car::leveneTest(Taxonomic_evenness ~ factor(Biomass_group), data = SEM_data[!is.na(SEM_data$Biomass_group), ])
t_test_result_evenness <- t.test(Taxonomic_evenness ~ Biomass_group, data = SEM_data, var.equal = FALSE, na.action = na.omit)
t_test_result_evenness

round(mean(SEM_data$Taxonomic_evenness[SEM_data$Biomass_group == "≥75th percentile"], na.rm = TRUE), 2)
round(sd(SEM_data$Taxonomic_evenness[SEM_data$Biomass_group == "≥75th percentile"], na.rm = TRUE), 2)
round(mean(SEM_data$Taxonomic_evenness[SEM_data$Biomass_group == "≤25th percentile"], na.rm = TRUE), 2)
round(sd(SEM_data$Taxonomic_evenness[SEM_data$Biomass_group == "≤25th percentile"], na.rm = TRUE), 2)

#
round(mean(SEM_data$Taxonomic_richness[SEM_data$Biomass_group == "≥75th percentile"], na.rm = TRUE), 2)
round(sd(SEM_data$Taxonomic_richness[SEM_data$Biomass_group == "≥75th percentile"], na.rm = TRUE), 2)
round(mean(SEM_data$Taxonomic_richness[SEM_data$Biomass_group == "≤25th percentile"], na.rm = TRUE), 2)
round(sd(SEM_data$Taxonomic_richness[SEM_data$Biomass_group == "≤25th percentile"], na.rm = TRUE), 2)

t_test_result_richness <- t.test(Taxonomic_richness ~ Biomass_group, data = SEM_data, var.equal = FALSE, na.action = na.omit)
t_test_result_richness

# Spatial pattern of evenness in the Adriatic Sea
cor.test(SEM_data_Adriatic$Taxonomic_evenness, SEM_data_Adriatic$Latitude, method = "spearman")
t_test_result_evenness <- t.test(Taxonomic_evenness ~ Latitude_group, data = SEM_data_Adriatic, var.equal = FALSE)
t_test_result_evenness

# Correlation between functional richness and Log biomass
GSA_Adriatic <- SEM_data %>% filter(GSA == "Adriatic Sea")
GSA_GOL <- SEM_data %>% filter(GSA == "Gulf of Lion")
GSA_Aegean <- SEM_data %>% filter(GSA == "Aegean Sea")
GSA_NSpain <- SEM_data %>% filter(GSA == "Northern Spain")

cor.test(GSA_Adriatic$Log_biomass, GSA_Adriatic$Functional_richness, method = "spearman")
cor.test(GSA_GOL$Log_biomass, GSA_GOL$Functional_richness, method = "spearman")
cor.test(GSA_Aegean$Log_biomass, GSA_Aegean$Functional_richness, method = "spearman")
cor.test(GSA_NSpain$Log_biomass, GSA_NSpain$Functional_richness, method = "spearman")

cor.test(SEM_data$Log_biomass, SEM_data$Functional_richness,  method = "spearman")

# plot relationships between Log biomass and functional richness per GSA
ggplot(SEM_data, aes(Log_biomass, Functional_richness)) +
  geom_point(aes(fill = GSA), shape = 21) +
  # geom_smooth per GSA
  geom_smooth(aes(color = GSA), method = "lm", se = FALSE) +
  theme_classic()

# plot relationships between functional richness and taxonomic richness per GSA
ggplot(SEM_data, aes(Taxonomic_richness, Functional_richness)) +
  geom_point(aes(fill = GSA), shape = 21) +
  # geom_smooth per GSA
  geom_smooth(aes(color = GSA), method = "lm", se = FALSE) +
  theme_classic()

cor.test(SEM_data$Taxonomic_richness, SEM_data$Functional_richness, method = "spearman")
cor.test(GSA_Adriatic$Taxonomic_richness, GSA_Adriatic$Functional_richness, method = "spearman")

# Spatial pattern of functional richness in the Adriatic Sea
cor.test(SEM_data_Adriatic$Functional_richness, SEM_data_Adriatic$Latitude, method = "spearman")

# Correlation between functional evenness and Log biomass
cor.test(SEM_data$Functional_evenness, SEM_data$Log_biomass, method = "spearman")

# Correlation between CWM length and Log biomass
cor.test(SEM_data$CWM_length, SEM_data$Log_biomass, method = "spearman")

# Comparing CWM length from hauls in the top biomass quartile (≥75th percentile) vs. the bottom biomass quartile (≤25th percentile)
car::leveneTest(CWM_length ~ factor(Biomass_group), data = SEM_data[!is.na(SEM_data$Biomass_group), ])
t_test_result_CWM_length <- t.test(CWM_length ~ factor(Biomass_group), data = SEM_data[!is.na(SEM_data$Biomass_group), ], var.equal = TRUE)

mean(SEM_data$CWM_length[SEM_data$Biomass_group == "≥75th percentile"], na.rm = TRUE)
mean(SEM_data$CWM_length[SEM_data$Biomass_group == "≤25th percentile"], na.rm = TRUE)
sd(SEM_data$CWM_length[SEM_data$Biomass_group == "≥75th percentile"], na.rm = TRUE)
sd(SEM_data$CWM_length[SEM_data$Biomass_group == "≤25th percentile"], na.rm = TRUE)

# Spatial pattern of CWM length in the Adriatic Sea
cor.test(SEM_data_Adriatic$CWM_length, SEM_data_Adriatic$Latitude, method = "spearman")

# Student's t-test to compare CWM length ≥43°N in Adriatic Sea vs. <43°N in Adriatic Sea
car::leveneTest(CWM_length ~ Latitude_group, data = SEM_data_Adriatic)
t_test_result_CWM_length_latitude <- t.test(CWM_length ~ Latitude_group, data = SEM_data_Adriatic, var.equal = TRUE)

mean(SEM_data_Adriatic$CWM_length[SEM_data_Adriatic$Latitude_group == "≥43°N"], na.rm = TRUE)
mean(SEM_data_Adriatic$CWM_length[SEM_data_Adriatic$Latitude_group == "<43°N"], na.rm = TRUE)
sd(SEM_data_Adriatic$CWM_length[SEM_data_Adriatic$Latitude_group == "≥43°N"], na.rm = TRUE)
sd(SEM_data_Adriatic$CWM_length[SEM_data_Adriatic$Latitude_group == "<43°N"], na.rm = TRUE)

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
View(metrics_summary)

# Compute statistical differences between GSAs for each metric using ANOVA and Tukey's HSD test
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

ggsave("Figures/Figure_S2_gsa_metrics.png", plot = combined_plot, width = 12, height = 17, dpi = 300)
# Figure caption: "Figure S2. Mean values of log biomass, taxonomic richness, taxonomic evenness, functional richness, functional evenness, and log CWM length across the four GSAs.
# Error bars represent standard deviations.
# Different letters indicate significant differences among GSAs based on Games–Howell post hoc pairwise comparisons (p < 0.05)."