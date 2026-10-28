# ==============================================================================
# Summary statistics and among-GSA comparisons
#
#   Correlations, t-tests, Levene tests and Welch's ANOVA with Games-Howell comparisons of
#   the biodiversity and biomass metrics.
#
# ==============================================================================

# ==============================================================================
# 1. Load libraries and prepare the environment
# ==============================================================================

# Load libraries
source("./0-Load libraries.R")

rm(list = ls())

load("Data/SEM_data_all_GSA.RData")
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

round(mean(SEM_data$Taxonomic_evenness[SEM_data$Biomass_group == "≥75th percentile"], na.rm = TRUE), 2)
round(sd(SEM_data$Taxonomic_evenness[SEM_data$Biomass_group == "≥75th percentile"], na.rm = TRUE), 2)
round(mean(SEM_data$Taxonomic_evenness[SEM_data$Biomass_group == "≤25th percentile"], na.rm = TRUE), 2)
round(sd(SEM_data$Taxonomic_evenness[SEM_data$Biomass_group == "≤25th percentile"], na.rm = TRUE), 2)

round(mean(SEM_data$Taxonomic_richness[SEM_data$Biomass_group == "≥75th percentile"], na.rm = TRUE), 2)
round(sd(SEM_data$Taxonomic_richness[SEM_data$Biomass_group == "≥75th percentile"], na.rm = TRUE), 2)
round(mean(SEM_data$Taxonomic_richness[SEM_data$Biomass_group == "≤25th percentile"], na.rm = TRUE), 2)
round(sd(SEM_data$Taxonomic_richness[SEM_data$Biomass_group == "≤25th percentile"], na.rm = TRUE), 2)

t_test_result_richness <- t.test(Taxonomic_richness ~ Biomass_group, data = SEM_data, var.equal = FALSE, na.action = na.omit)

# Spatial pattern of evenness in the Adriatic Sea
cor.test(SEM_data_Adriatic$Taxonomic_evenness, SEM_data_Adriatic$Latitude, method = "spearman")
t_test_result_evenness <- t.test(Taxonomic_evenness ~ Latitude_group, data = SEM_data_Adriatic, var.equal = FALSE)

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
  geom_smooth(aes(color = GSA), method = "lm", se = FALSE) +
  theme_classic()

# plot relationships between functional richness and taxonomic richness per GSA
ggplot(SEM_data, aes(Taxonomic_richness, Functional_richness)) +
  geom_point(aes(fill = GSA), shape = 21) +
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