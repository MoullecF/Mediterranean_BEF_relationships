# ==============================================================================
# Build the SEM dataset
# ==============================================================================

# ==============================================================================
# 1. Load libraries and prepare the environment
# ==============================================================================

# Load libraries
source(here::here("R/0-Load libraries.R"))

rm(list = ls())

load(here::here("Tools/make_correlation_plot_function.RData"))

# ==============================================================================
# 2. Load and merge data
# ==============================================================================

# GSAs to process
GSA <- c(6, 7, 17, 18, 22, 23)
load(paste0(here::here("Data"), "/community_GSA_", paste0(GSA, collapse = "_"), ".RData"))
load(paste0(here::here("Outputs"), "/diversity_GSA_", paste0(GSA, collapse = "_"), ".RData"))
load(paste0(here::here("Data"), "/Environment/environmental_data_GSA_", paste0(GSA, collapse = "_"), ".RData"))

haul_biomass <- community %>%
    summarise(wgt_cpua = sum(wgt_cpua, na.rm = TRUE), .by = haul_id)

haul_meta <- community %>%
    distinct(haul_id, month, year, longitude, latitude) %>%
    mutate(GSA = as.numeric(stringr::str_extract(haul_id, "^[0-9]+")))

data <- diversity %>%
    left_join(haul_biomass, by = "haul_id") %>%
    left_join(haul_meta, by = "haul_id") %>%
    left_join(environmental_data, by = "haul_id") %>%
    filter(month %in% c(5, 6, 7, 8)) %>%
    mutate(
      Biomass = wgt_cpua,
      log_biomass = log(Biomass),
      GSA = GSA) %>%
    dplyr::select(haul_id, tric, teve, fric, feve, Length_max_cwm, Biomass, log_biomass, 
    depth, sst_mean, sst_sd, sbt_mean, sbt_sd, so_mean, so_sd, chla_mean, FPI_OTB, FPI_PS, FPI_SSF, gravity, year, GSA, longitude, latitude) %>%
    as.data.frame()

# Check NA values in the data frame
na_summary <- sapply(data, function(x) sum(is.na(x)))

# Plot NA on a map to see if they are spatially structured
na_data <- data %>% filter(if_any(everything(), is.na))
ggplot(na_data, aes(x = longitude, y = latitude)) +
  borders("world", colour = "gray50", fill = "gray80") +
  coord_cartesian(xlim = c(-6, 36), ylim = c(30, 46)) +
  geom_point(color = "red", size = 0.1) +
  theme_bw() +
  labs(title = "Spatial distribution of NA values", x = "Longitude", y = "Latitude")
# NAs are in the Gulf of Euboea (no env  data) so ok

# Remove hauls with NA values in any of the variables of interest and create a key for coordinates to check for duplicates
data <- data %>% 
  drop_na() %>%
  mutate(coords_key = paste(longitude, latitude, sep = "_"))

# Check for duplicates in data coordinates
if (any(duplicated(data$coords_key))) {
  message(paste0("There are ", sum(duplicated(data$coords_key)), " duplicate coordinates."))
  message("[Spatial] Duplicate coordinates detected — applying small jitter (1e-5 degrees).")

# Add very small jitter to duplicate coordinates
  data <- data %>%
    group_by(coords_key) %>%
    mutate(
      longitude = if_else(duplicated(coords_key) | duplicated(coords_key, fromLast = TRUE),
                          longitude + runif(n(), -1e-5, 1e-5),
                          longitude),
      latitude = if_else(duplicated(coords_key) | duplicated(coords_key, fromLast = TRUE),
                         latitude + runif(n(), -1e-5, 1e-5),
                         latitude)
    ) %>%
    ungroup()

  # Check again for duplicates
  data <- data %>%
    mutate(coords_key = paste(longitude, latitude, sep = "_"))
  if (any(duplicated(data$coords_key))) {
    message("There are still duplicate coordinates after jittering !")
    message(paste0("There are still ", sum(duplicated(data$coords_key)), " duplicate coordinates after jittering."))
  } else {
    message("[Spatial] No more duplicate coordinates after jittering.")
  }
}

# Identify outliers for low and high values of Biomass (total biomass)
options(max.print = 2000)
print(quantile(data$Biomass, probs = seq(0, 1, 0.0005))) 

# Remove biomass outliers, log-transform chla_mean and remove hauls with CWM of max length > 250 cm (as only one rare species is driving the index (HEXAGRI))
data <- data %>%
  filter(Biomass > 15 & Biomass < 5700) %>% # removed 0.5% of lowest and highest values (1% in total)
  mutate(log_chla_mean = log(chla_mean)) %>%
  filter(Length_max_cwm <= 250)

# Number of hauls per GSA after filtering, group 17 and 18 and 22 and 23 together as they are the "same" sub-regions
data %>%
  mutate(GSA_grouped = dplyr::recode(GSA,
    `6` = "6",
    `7` = "7",
    `17` = "17_18",
    `18` = "17_18",
    `22` = "22_23",
    `23` = "22_23"
  )) %>%
  group_by(GSA_grouped) %>%
  summarise(n_hauls = n())

# ==============================================================================
# 3. Collinearity screening
# ==============================================================================

### Test for correlation between covariates
df_cor1 = data[, c("tric", "teve", "feve", "fric",
                  "Length_max_cwm", "log_biomass", "depth",
                  "sst_mean", "sst_sd", "sbt_mean", "sbt_sd", "so_mean", 
                  "so_sd", "chla_mean", "FPI_OTB", "FPI_PS", "FPI_SSF", "gravity")]

cor1 <- make_correlation_plot(df_cor = df_cor1, threshold = 0.7, method = "pearson") 
usdm::vifstep(as.data.frame(df_cor1), th=5, keep= c("log_biomass"), size = 10100) # keep log_biomass as it is the response variable

#  sbt_sd and FPI_PS removed
df_cor2 = data[, c("tric", "teve", "feve", "fric", "Length_max_cwm", "log_biomass",
                    "depth", "sst_mean", "sst_sd", "sbt_mean", "so_mean", 
                  "so_sd", "chla_mean", "FPI_OTB", "FPI_SSF", "gravity")]


corr_fig <- ggcorrplot(cor(df_cor2), type = "lower", lab = TRUE, method="circle", colors = c("red", "white", "blue"), digits = 2, legend.title = "Pearson\nCorrelation", ggtheme=theme_bw())
usdm::vifstep(as.data.frame(df_cor2), th=5, keep= c("log_biomass"), size = 10100) # keep log_biomass as it is the response variable
# No variable from the 16 input variables has collinearity problem. 

# ==============================================================================
# 4. Build and save the SEM dataset
# ==============================================================================

# Create data frame for SEM with only the variables of interest
SEM_data <- data %>%
    dplyr::select(c("log_biomass", # biomass
                    "tric", "teve", # taxonomic diversity metrics
                    "feve", "fric", # functional diversity metrics
                    "Length_max_cwm", # community-weighted mean index
                    "depth", "sst_mean", "sst_sd", "sbt_mean", "so_mean", "so_sd", "chla_mean", # environmental variables
                    "FPI_OTB", "FPI_SSF", "gravity", # fishing pressure and human impacts
                    "year", "GSA","longitude", "latitude" # for random effects and spatial autocorrelation
)) %>%
# Change names of GSA to be more readable
  dplyr::mutate(GSA = case_when(
    GSA == 6 ~ "Catalan_Sea",
    GSA == 7 ~ "Gulf_of_Lion",
    GSA == 17 ~ "Adriatic_Sea",
    GSA == 18 ~ "Adriatic_Sea",
    GSA == 22 ~ "Aegean_Sea",
    GSA == 23 ~ "Aegean_Sea")) %>%
dplyr::mutate(GSA = as.factor(GSA)) %>%
as.data.frame()

# Save the data frame for SEM
# save(SEM_data, file = here::here("Data/SEM_data_all_GSA.RData"))
