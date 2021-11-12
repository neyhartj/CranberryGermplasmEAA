# Germplasm environmental association
#
# Download bioclim data for environmental associations
#

# Load packages
library(raster)
library(tidyverse)
library(readxl)


# Set directories
data_dir <- here::here("data")
results_dir <- here::here("results")

# External directories
# CranberryLab directory
path_split <- str_split(string = here::here(), pattern = "/")[[1]]
cran_dir <- paste0(path_split[seq_len(str_which(string = path_split, pattern = "CranberryLab"))], collapse = "/")
# Directory of metadata
meta_dir <- file.path(cran_dir, "Breeding/prior2021/Populations/GermplasmCollection")
# Directory for environmental data
env_dir <- file.path(cran_dir, "EnvironmentalData")


# Resolution at which to aggregate soil data (in minutes)
soil_data_aggr_res <- 2.5



# Read in data ------------------------------------------------------------

# Read in population metadata - local
load(file.path(data_dir, "germplasm_metadata.RData"))
# Load the local marker data
load(file.path(data_dir, "gc_marker_data.RData"))


# Create abbreviations for the origin names
location_metadata <- germ_meta %>%
  # filter for missing origin or coordinates; select only genotypes in K
  filter(!is.na(origin_name),
         !is.na(origin_latitude),
         individual %in% row.names(K_wild)) %>%
  distinct(origin_name, origin_latitude, origin_longitude) %>%
  mutate(origin_name1 = str_remove(origin_name, "west branch, ")) %>%
  separate(origin_name1, c("location", "state"), sep = ", ") %>%
  mutate(origin_abbr = abbreviate(location, 3),
         origin_abbr = paste0(toupper(origin_abbr), state)) %>%
  select(location, origin_name, state, location_abbr = origin_abbr,
         latitude = origin_latitude, longitude = origin_longitude)

# Set the min max of lat/long for subsetting data
lat_range <- range(location_metadata$latitude, na.rm = TRUE)
long_range <- range(location_metadata$longitude, na.rm = TRUE)



# Get altitude and WorldClim data ------------------------------------------------------

# Many lines of this code are from here:
# https://github.com/MorrellLAB/Env_Assoc/blob/master/script/GWAS/worldclim.miner3.r

# Get worldclim data, bioclimatic variables at 5 minute resolution (~9 sq km)
#
#
# # This is only to be used once; Use the next command to read in data
# worldclim_dat <- getData(name = "worldclim", var = "bio", res = "5",
#                          path = data_dir)
#
# download.file(url = "https://biogeo.ucdavis.edu/data/diva/msk_alt/USA_msk_alt.zip",
#               destfile = file.path(env_dir, "Elevation/USA_msk_alt.zip"))
# download.file(url = "https://biogeo.ucdavis.edu/data/diva/msk_alt/CAN_msk_alt.zip",
#               destfile = file.path(env_dir, "Elevation/CAN_msk_alt.zip"))

# Read in the elevation data
usa_alt_data <- getData(name = "alt", country = "USA", download = FALSE, keepzip = TRUE,
                        path = file.path(env_dir, "Elevation"))
can_alt_data <- getData(name = "alt", country = "CAN", path = file.path(env_dir, "Elevation"),
                        keepzip = TRUE)


# Combine; usa[[1]] is continental
alt_data <- raster::merge(x = usa_alt_data[[1]], y = can_alt_data)


# Read in already downloaded worldclim data
worldclim_dat <- getData(name = "worldclim", download = FALSE, var = "bio", res = "5",
                         path = file.path(env_dir, "WorldClim"))

# Crop this data
long_range1 <- range(pretty(long_range))
lat_range1 <- range(pretty(lat_range))
extent1 <- extent(list(x = long_range1, y = lat_range1))
worldclim_dat1 <- crop(x = worldclim_dat, y = extent1)
alt_data1 <- crop(x = alt_data, y = extent1)

# Names of bioclim variables
environ_vars <- c(
  elevation = "Elevation",
  latitude = "Latitude",
  longitude = "Longitude",
  bio1 = 'Annual Mean Temperature',
  bio2 = 'Mean Diurnal Range (Mean of monthly (max temp - min temp))',
  bio3 = 'Isothermality (bio2/bio7) (* 100)',
  bio4 = 'Temperature Seasonality (standard deviation *100)',
  bio5 = 'Max Temperature of Warmest Month',
  bio6 = 'Min Temperature of Coldest Month',
  bio7 = 'Temperature Annual Range (bio5-bio6)',
  bio8 = 'Mean Temperature of Wettest Quarter',
  bio9 = 'Mean Temperature of Driest Quarter',
  bio10 = 'Mean Temperature of Warmest Quarter',
  bio11 = 'Mean Temperature of Coldest Quarter',
  bio12 = 'Annual Precipitation',
  bio13 = 'Precipitation of Wettest Month',
  bio14 = 'Precipitation of Driest Month',
  bio15 = 'Precipitation Seasonality (Coefficient of Variation)',
  bio16 = 'Precipitation of Wettest Quarter',
  bio17 = 'Precipitation of Driest Quarter',
  bio18 = 'Precipitation of Warmest Quarter',
  bio19 = 'Precipitation of Coldest Quarter'
)

# Combine into tidy df with classifiers for each variable
environ_vars_df <- tibble(variable = names(environ_vars), full_name = environ_vars) %>%
  mutate(class = case_when(variable %in% paste0("bio", 1:11) ~ "temperature",
                           variable %in% paste0("bio", 12:19) ~ "precipitation",
                           TRUE ~ "geography"))



# Subset data for each coordinate
location_metadata1 <- location_metadata %>%
  filter_at(vars(latitude, longitude), all_vars(!is.na(.))) %>%
  mutate(data = list(NULL))

# Iterate over rows
for (i in seq_len(nrow(location_metadata1))) {

  # Get latitude/longitude
  lat_i <- location_metadata1$latitude[i]
  long_i <- location_metadata1$longitude[i]

  # Build extent
  geo_add <- 1e-5
  extent_i <- extent(list(x = c(long_i - geo_add, long_i + geo_add),
                          y = c(lat_i - geo_add, lat_i + geo_add)))

  # Crop the data
  worldclim_dat_i <- crop(x = worldclim_dat1, y = extent_i)
  alt_data_i <- crop(x = alt_data1, y = extent_i)

  # Create an empty data.frame
  bioclim_dat <- t(as.matrix(worldclim_dat_i)) %>%
    as.data.frame() %>%
    rownames_to_column("bioc_var") %>%
    rename(value = V1) %>%
    as_tibble()

  # Add elevation
  alt_i <- (as.matrix(alt_data_i)[,]) %>% ifelse(is.na(.), 0, .)
  bioclim_dat1 <- bioclim_dat %>%
    add_row(bioc_var = c("elevation", "latitude", "longitude"),
            value = c(as.matrix(alt_data_i)[,], lat_i, long_i))

  # Add this tibble to the larger data.frame
  location_metadata1$data[[i]] <- bioclim_dat1

}

# Unnest
location_metadata2 <- unnest(location_metadata1, data) %>%
  spread(bioc_var, value) %>%
  # Convert temperature data
  mutate_at(all_of(subset(environ_vars_df, class == "temperature", variable, drop = TRUE)), ~ . / 10)

location_worldclim_data <- location_metadata2

# Save this
save("location_worldclim_data", "environ_vars_df", file = file.path(data_dir, "germplasm_origin_worldclim_data.RData"))








