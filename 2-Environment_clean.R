rm(list = ls())

## Library ----
load("tools/install_load_function.RData")

# Install/load required packages
install_load_function(
  
  # list of required packages
  c("dplyr",
    "tidyr",
    "ggplot2",
    "reticulate", # Python interface
    "copernicusR",
    "usethis", # edit R profile or environment
    "lubridate", # work with raster dates
    "sf",
    "terra",
    "gfwr"
  ))

# Set GSA number
GSA <- c(6, 7, 17, 18, 22, 23)

load(paste0("Data/community_GSA_", paste0(GSA, collapse = "_"), ".RData"))
load(paste0("Outputs/diversity_GSA_", paste0(GSA, collapse = "_"), ".RData"))

# Merge haul-level diversity with haul-level total biomass from community
haul_biomass <- community %>%
  summarise(wgt_cpua = sum(wgt_cpua, na.rm = TRUE), .by = haul_id)

haul_meta <- community %>%
  distinct(haul_id, month, year, depth, longitude, latitude)

# total Number of unique hauls and number of hauls with tric < 5
total_hauls <- n_distinct(diversity$haul_id) # 11093 unique hauls
hauls_tric_lt_5 <- sum(diversity$tric < 5, na.rm = TRUE) # 0!

data <- diversity %>%
  left_join(haul_biomass, by = "haul_id") %>%
  left_join(haul_meta, by = "haul_id") %>%
  filter(month %in% c(5, 6, 7, 8)) %>% # select only month 5 to 8. Remove 928 hauls
  filter(tric > 4) %>% # remove hauls with tric < 5 (no effect on the dataset)
  as.data.frame()

### Download Copernicus data ----

# Install miniconda (Python) -
# https://www.anaconda.com/docs/getting-started/miniconda/main
# install_miniconda()

# Check that miniconda is installed
#conda_list()

# Then install the copernicus marine toolbox
# reticulate::py_install("copernicusmarine")

# Tell R which Python environment to use
use_condaenv("r-reticulate", required = TRUE)
py_config() # check path and version

# Access the data using your credentials (Copernicus account)
# Register at: https://data.marine.copernicus.eu/register?pk_vid=11ceae9d3160528d17743415560fb363
setup_copernicus(username = "fabien.moullec@umontpellier.fr",
                 password = "C7FR!b75z!zrn.A",
                 store_credentials = TRUE)

# Define spatial and temporal extent of the data to download
ye <- range(data$year)
la <- range(data$latitude)
lo <- range(data$longitude)

# Download temperature data (surface and bottom) from MEDSEA_MULTIYEAR_PHY_006_004
f_1 <- copernicus_download(
  dataset_id = "cmems_mod_med_phy-temp_my_4.2km_P1M-m", # monthly
  variables = c("thetao", "bottomT"), # sst, sbt
  bbox = c(lo[1] - 0.3, 
           lo[2] + 0.3,
           la[1] - 0.3,
           la[2] + 0.3),
  start_date = paste0(ye[1] - 1, "-01-01"),
  end_date = paste0(ye[2], "-12-31"),
  depth = c(1.02 , 1.02 ), # first depth layer
  dataset_version = "202511")

# Move to data folder
file.rename(f_1, paste0("Data/Environment/temperature_data_GSA_",paste0(GSA, collapse = "_"),".nc"))

# Read data as raster 
raster_t_data <- terra::rast(paste0("Data/Environment/temperature_data_GSA_",paste0(GSA, collapse = "_"),".nc"))

# Check layers
names(raster_t_data)

# Visual check
terra::plot(raster_t_data[[1:12]])

# Download salinity data from MEDSEA_MULTIYEAR_PHY_006_004
f_2 <- copernicus_download(dataset_id = "cmems_mod_med_phy-sal_my_4.2km_P1M-m", # monthly
                           variables = c("so"), # salinity
                           bbox = c(lo[1] - 0.3, 
                                    lo[2] + 0.3,
                                    la[1] - 0.3,
                                    la[2] + 0.3),
                           start_date = paste0(ye[1] - 1, "-01-01"),
                           end_date = paste0(ye[2], "-12-31"),
                           depth = c(1.02 , 1.02 ), # first depth layer
                           dataset_version = "202511")


# Move to data folder
file.rename(f_2, paste0("Data/Environment/salinity_data_GSA_",paste0(GSA, collapse = "_"),".nc"))

# Read data as raster 
raster_s_data <- terra::rast(paste0("Data/Environment/salinity_data_GSA_",paste0(GSA, collapse = "_"),".nc"))

# Check layers
names(raster_s_data)

# Visual check
terra::plot(raster_s_data[[1:12]])

# Download biogeochemical data from Global Ocean Colour
# https://doi.org/10.48670/moi-00281
f_3 <- copernicus_download(
  dataset_id = "cmems_obs-oc_glo_bgc-plankton_my_l4-multi-4km_P1M", # monthly
  variables = c("CHL"), # Chla [mg/m3]
  bbox = c(lo[1] - 0.3, 
           lo[2] + 0.3,
           la[1] - 0.3,
           la[2] + 0.3),
  start_date = paste0(ye[1] - 1, "-01-01"),
  end_date = paste0(ye[2], "-12-31"),
  dataset_version = "202603" # until end of 2023
)

# Move to data folder (may fail across drives)
file.rename(f_3, paste0("Data/Environment/chlorophyll_data_GSA_",paste0(GSA, collapse = "_"),".nc"))

# Read data as raster
chla_stack <- rast(paste0("Data/Environment/chlorophyll_data_GSA_",paste0(GSA, collapse = "_"),".nc"))

# Check layers
names(chla_stack)

# Visual check
terra::plot(chla_stack[[1:12]])

## Download Global Fishing Watch data ----

# Install gfwr package
#remotes::install_github("GlobalFishingWatch/gfwr",
#                        dependencies = TRUE)

# To access the API, a token is needed. It can be requested here:
# https://gateway.api.globalfishingwatch.org/v3/auth?client=gfw&callback=https%3A%2F%2Fglobalfishingwatch.org%2Four-apis%2Ftokens&locale=en

# Then, the token must be saved to the .Renviron file using 
#usethis::edit_r_environ() # write GFW_TOKEN="YOUR_TOKEN" and save .Renviron

# Restart R for changes to take effect

# Check that the token variable has been set
#Sys.getenv("GFW_TOKEN") # should return the token

# Create the object key with the token
#key <- gfw_auth()

# Create an sf object to define the spatial extent of the data
#ye <- range(data$year)
#la <- range(data$latitude)
#lo <- range(data$longitude)

#polygon <- sf::st_bbox(c(xmin = lo[1] - 1,
  #                        xmax = lo[2] + 1,
  #                        ymin = la[1] - 1,
  #                        ymax = la[2] + 1),
  #                      crs = 4326) %>%
  # sf::st_as_sfc() %>%
  # sf::st_as_sf()

# Only one year of data can be imported at a time, so a loop is necessary
# years <- (ye[1] - 1):ye[2]

# fishing_list <- vector("list", length(years))

# for (i in seq_along(years)) {
  
#   cat("Downloading fishing data for year", years[i], "...\n")
#   yr <- years[i]
  
#   fishing_list[[i]] <- gfw_ais_fishing_hours(
#     spatial_resolution = "HIGH",    # 0.01 degree or ~ 1km
#     temporal_resolution = "MONTHLY",
#     group_by = "GEARTYPE",
#     start_date = paste0(yr, "-01-01"),
#     end_date = paste0(yr, "-12-31"),
#     region_source = "USER_SHAPEFILE",
#     region = polygon
#   )
# }

# # Bind all data frames
# fishing_data <- do.call(rbind, fishing_list)

# # Save
# save(fishing_data, file = paste0("Data/Environment/fishing_data_GSA_",GSA,".RData"))

# # What gears are used in the area?
# prop.table(table(fishing_data$geartype)) # no GEAR TYPE !

# # Keep only trawling and prepare for raster format
# trawling_data <- fishing_data %>%
#   #filter(geartype == "trawlers") %>% # focus on trawling (60% of observations)
#   mutate(time = as.Date(paste0(`Time Range`, "-01"))) %>% # prepare time column
#   group_by(time, Lat, Lon) %>%
#   summarize(fishing_hours = sum(`Apparent Fishing Hours`, na.rm = TRUE),
#             .groups = "drop") 

# # Prepare a raster grid to extract values
# r_grid <- rast(
#   xmin = min(trawling_data$Lon),
#   xmax = max(trawling_data$Lon),
#   ymin = min(trawling_data$Lat),
#   ymax = max(trawling_data$Lat),
#   resolution = 0.01  # fishing_data resolution
# )

# # Split data by time (month and year)
# time_list <- split(trawling_data, trawling_data$time) # these will be the layers

# # Rasterise each time layer
# r_list <- lapply(time_list, function(df) {
  
#   pts <- vect(df, geom = c("Lon", "Lat"), crs = "EPSG:4326")
  
#   rasterize(
#     pts,
#     r_grid,
#     field = "fishing_hours",
#     fun = "sum"   # sum if multiple points per cell
#   )
# })

# # Stack the layers and name them
# fishing_stack <- rast(r_list)

# time(fishing_stack) <- as.Date(names(time_list))

# # Check
# fishing_stack

# terra::plot(log(fishing_stack[[100]]), main = paste0("Fishing hours in ", names(fishing_stack)[100], " (log scale)")) # check a layer

# # Save intermediate
# writeRaster(fishing_stack, file = paste0("Data/Environment/fish_stack_GSA_",GSA,".tiff"), overwrite = TRUE)


## Compute mean, sd, and sum values ----

# Load rasters 
temperature_stack <- terra::rast(paste0("Data/Environment/temperature_data_GSA_",paste0(GSA, collapse = "_"),".nc"))
chla_stack <- terra::rast(paste0("Data/Environment/chlorophyll_data_GSA_",paste0(GSA, collapse = "_"),".nc"))
salinity_stack <- terra::rast(paste0("Data/Environment/salinity_data_GSA_",paste0(GSA, collapse = "_"),".nc"))
FPI_OTB <- get(load("Data/Environment/FPI_OTB_stack.RData")) %>% terra::rast(FPI_OTB)
FPI_PS <- get(load("Data/Environment/FPI_PS_stack.RData")) %>% terra::rast(FPI_PS)
FPI_SSF <- get(load("Data/Environment/FPI_SSF_stack.RData")) %>% terra::rast(FPI_SSF)
Gravity_layer <- get(load("./Data/Environment/r_gravity_Med.RData")) %>% terra::rast(Gravity_layer)

# Subset physical data to compute each variable on its own
sst_stack <- temperature_stack[[grep("thetao", names(temperature_stack))]] # sst
sbt_stack <- temperature_stack[[grep("bottomT", names(temperature_stack))]]# sbt
so_stack <- salinity_stack[[grep("so", names(salinity_stack))]] # salinity
Gravity_layer <- log(Gravity_layer) # log-transform gravity

# Output data frame
environmental_data <- data %>%
  group_by(haul_id, year, latitude, longitude) %>%
  summarize(depth = median(depth), .groups = "drop")

# CRS for extraction (equal area)
#y_center <- round(mean(range(environmental_data$latitude, na.rm = TRUE)))
#x_center <- round(mean(range(environmental_data$longitude, na.rm = TRUE)))

# equal_crs <- paste0("+proj=laea +lat_0=",
#                     y_center,
#                     " +lon_0=",
#                     x_center,
#                     " +datum=WGS84 +units=km +no_defs")

# Haul point positions for extraction
haul_points <- environmental_data %>% 
  vect(geom = c("longitude", "latitude"),
       crs = "EPSG:4326") #%>% # should be the correct crs
  #project(equal_crs)

# Month to divide years
earliest_sampling_month <- min(as.numeric(unique(data$month)))
  
# Loop to compute yearly mean and sd, project, and extract
stacks <- list(
  sst  = "sst_stack",
  sbt  = "sbt_stack",
  so   = "so_stack",
  chla = "chla_stack"
)

for (v in names(stacks)) {
  
  r <- get(stacks[[v]])
  
  # Shift year to identify the 12 months before sampling ----
  dates <- time(r)
  yrs   <- as.integer(format(dates, "%Y"))
  mons  <- as.integer(format(dates, "%m"))
  
  yrs_new <- ifelse(mons >= earliest_sampling_month,
                    yrs + 1, yrs)
  
  # Compute yearly mean and sd ----
  r_mean <- tapp(r, index = yrs_new, fun = mean, na.rm = TRUE)
  r_sd   <- tapp(r, index = yrs_new, fun = sd,   na.rm = TRUE)
  
  # Assign time ----
  time(r_mean) <- as.numeric(gsub("X", "", names(r_mean)))
  time(r_sd)   <- as.numeric(gsub("X", "", names(r_sd)))
  
  # Extract mean value within a 2 km buffer of each haul point ----
  buf_2km <- terra::buffer(haul_points, width = 2000)
  mean_vals <- terra::extract(r_mean, buf_2km, fun = mean, na.rm = TRUE)[, -1]
  sd_vals   <- terra::extract(r_sd,   buf_2km, fun = mean, na.rm = TRUE)[, -1]
  
  # Fill remaining NAs using the mean of non-NA cells within 10 km ----
  # (2 km buffer had no non-NA cells at these hauls, e.g. near coastlines;
  # widen the search radius to 10 km)
  for (j in seq_len(ncol(mean_vals))) {
  
  rows_j <- which(is.na(mean_vals[, j]))
  if (length(rows_j) == 0) next
  
  buf_j <- terra::buffer(haul_points[rows_j], width = 10000)
  fill_vals <- terra::extract(r_mean[[j]], buf_j, fun = mean, na.rm = TRUE)
  
  mean_vals[rows_j, j] <- fill_vals[, 2]
}
  
  # Fill remaining NAs using the mean of non-NA cells within 10 km ----
  for (j in seq_len(ncol(sd_vals))) {
  
  rows_j <- which(is.na(sd_vals[, j]))
  if (length(rows_j) == 0) next
  
  buf_j <- terra::buffer(haul_points[rows_j], width = 10000)
  fill_vals <- terra::extract(r_sd[[j]], buf_j, fun = mean, na.rm = TRUE)
  
  sd_vals[rows_j, j] <- fill_vals[, 2]
}
  
  # Match raster layers to haul year ----
  layer_years <- as.numeric(gsub("X", "", names(r_mean)))
  idx <- match(environmental_data$year, layer_years)
  
  # Add new variables to data ----
  environmental_data[[paste0(v, "_mean")]] <- 
    mean_vals[cbind(seq_len(nrow(mean_vals)), idx)]
  
  environmental_data[[paste0(v, "_sd")]] <- 
    sd_vals[cbind(seq_len(nrow(sd_vals)), idx)]
}

# Check if there are any NAs in the extracted values
colSums(is.na(environmental_data))

load("Tools/plot_env_na.RData")
plot_env_na(environmental_data, raster = sst_stack[[1]], raster_name = "SST mean (°C)", env_vars = c("sst_mean"), point_size = 1.8, point_alpha = 1, buf = 1)

# Extract FPI values for each gear type for the year before
stacks <- list(FPI_OTB = "FPI_OTB", FPI_PS = "FPI_PS", FPI_SSF = "FPI_SSF")

environmental_data$FPI_OTB <- NA
environmental_data$FPI_PS <- NA
environmental_data$FPI_SSF <- NA

for (i in 1:dim(environmental_data)[1]) {
  
  cat("Extracting FPI values for haul", i, "...\n")
  
  for (v in names(stacks)) {
    
    r <- get(stacks[[v]])
    
    # Point and buffers at haul location
    point <- vect(environmental_data[i, c("longitude", "latitude")], 
                  geom = c("longitude", "latitude"), 
                  crs = "EPSG:4326")
    year <- environmental_data$year[i]

    r_year <- r[[grep(paste0(year - 1), names(r))]] # select layer for the year before

    # Extract mean value within 2 km first
    buf_2km <- terra::buffer(point, width = 2000)
    val <- terra::extract(r_year, buf_2km, fun = mean, na.rm = TRUE)[, 2]

    # Fill remaining NAs using the mean of non-NA cells within 10 km
    if (is.na(val)) {
      buf_10km <- terra::buffer(point, width = 10000)
      val <- terra::extract(r_year, buf_10km, fun = mean, na.rm = TRUE)[, 2]
    }

    # Transform to 0 if still NA (no fishing found even within 10 km)
    if (is.na(val)) {
      val <- 0
    }

    # Add to data frame
    environmental_data[i,stacks[[v]]] <- val
  }
}

# Add gravity index
gravity_points <- environmental_data %>%
  terra::vect(geom = c("longitude", "latitude"), crs = "EPSG:4326")

# Extract mean gravity value within a 2 km buffer of each haul point
gravity_buf_2km <- terra::buffer(gravity_points, width = 2000)
gravity_vals <- terra::extract(Gravity_layer, gravity_buf_2km, fun = mean, na.rm = TRUE)[, 2]

# Fallback: mean of non-NA cells within 10 km
na_idx <- is.na(gravity_vals)
if (any(na_idx)) {
  buf_na <- terra::buffer(gravity_points[na_idx], width = 10000)
  gravity_vals[na_idx] <- terra::extract(
    Gravity_layer,
    buf_na,
    fun    = mean,
    na.rm  = TRUE
  )[, 2]
}

environmental_data$gravity <- gravity_vals

plot_env_na(environmental_data, raster = Gravity_layer[[1]], raster_name = "Gravity index", env_vars = c("gravity"), point_size = 1.8, point_alpha = 1, buf = 1)

# Check if there are any NAs in the extracted values (should not be, but just in case)
colSums(is.na(environmental_data))

## Clean and save ----

# Final cleanup 
environmental_data <- environmental_data %>%
  dplyr::select(-c(year,
            latitude, # keep latitude and longitude for visual checks
            longitude, 
            chla_sd)) %>% # chla_sd makes little sense
  mutate(
    #across(c(sum_fishing), # NA fishing are actually 0
    #       ~ replace_na(., 0)), 
    log_chla_mean = log(chla_mean),# check your chla distribution
  log_depth = log(depth)) 

# Save
save(environmental_data, file = paste0("Data/Environment/environmental_data_GSA_", paste0(GSA, collapse = "_"),".RData"))