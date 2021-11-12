# Germplasm environmental association
#
# Combine the worldclim and soilgrids data into a single dataset
#

# Load packages
library(raster)
library(tidyverse)
library(readxl)
library(sf)


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



# Read in the acquired worldclim/soilgrid data ------------------------------------------------------------

# Load the worldclim data
load(file.path(data_dir, "germplasm_origin_worldclim_data.RData"))
load(file.path(data_dir, "germplasm_origin_soil_data.RData"))

# Spread out the aggregated soil data
location_soil_data_aggr1 <- location_soil_data_aggr %>%
  gather(variable, value, soil_vars_df$variable) %>%
  unite(var, variable, soil_layer, sep = "_") %>%
  spread(var, value)

# Combine the worldclim and soil data
location_bioclim_data <- full_join(location_worldclim_data, location_soil_data_aggr1)

# Edit the soil variable df
soil_vars_df1 <- soil_vars_df %>%
  crossing(., distinct(location_soil_data_aggr, soil_layer)) %>%
  mutate(variable = paste0(variable, "_", soil_layer),
         full_name = paste0(full_name, "_", soil_layer),
         full_name = str_to_title(str_replace_all(full_name, "_", " "))) %>%
  select(-soil_layer)

# Combine the list of variable names
bioclim_vars_df <- bind_rows(environ_vars_df, soil_vars_df1)

# Analyze latent factors -----------------------------------

# Center and scale
# Convert the data to a matrix
location_bioclim_data_scaled_mat <- location_bioclim_data %>%
  mutate_if(is.numeric, ~as.numeric(scale(.))) %>%
  select(location_abbr, starts_with("bio"), matches(paste(soil_vars_df$variable, collapse = "|"))) %>%
  as.data.frame() %>%
  column_to_rownames("location_abbr") %>%
  as.matrix()


## PCA first
pca_env = prcomp(location_bioclim_data_scaled_mat)
# Proportion of variation explained by each PC
pc_propvar <- (pca_env$sdev^2)/sum(pca_env$sdev^2)
# Cumulative sum
cumsum(pc_propvar)
# 81% by PCs 1-4
# 56% by PCs 1-2
# Just use 1-4

location_bioclim_data_scaled_mat_PCs = cbind(location_bioclim_data_scaled_mat, pca_env$x[,1:4])
round(cor(location_bioclim_data_scaled_mat_PCs)^2,digits=3)

## Only the first 3 PCs make sense

# How about ICA?
library(ica)

ic_env = icaimax(location_bioclim_data_scaled_mat, 3)
i_mat1 = ic_env$R
location_bioclim_data_scaled_mat_ICs = cbind(location_bioclim_data_scaled_mat, `colnames<-`(ic_env$S, paste0("IC", seq_along(ic_env$vafs))))
round(cor(location_bioclim_data_scaled_mat_ICs)^2,digits=3)

# Each variable has an IC that is decently correlated with it
ic_cor <- round(cor(location_bioclim_data_scaled_mat_ICs)^2,digits=3) %>%
  as.data.frame() %>%
  rownames_to_column("IC") %>%
  gather(variable, cor, -IC) %>%
  filter(str_detect(IC, "IC"), str_detect(variable, "bio|soil"))

# For each bio, find the IC most correlated
ic_cor %>%
  group_by(variable) %>%
  top_n(x = ., n = 1, wt = cor) %>%
  # Assign category
  left_join(., bioclim_vars_df) %>%
  # Split by IC
  split(.$IC)


# Add the ICs to the location_bioclim_data
eaa_environmental_data <- as.data.frame(location_bioclim_data_scaled_mat_ICs) %>%
  rownames_to_column("location_abbr") %>%
  select(location_abbr, starts_with("IC")) %>%
  left_join(location_bioclim_data, .)

eaa_environmental_vars <- bioclim_vars_df %>%
  add_row(variable = str_subset(colnames(location_bioclim_data_scaled_mat_ICs), "IC"),
          full_name = variable, class = "independent_component")

# Save this
save("eaa_environmental_data", "eaa_environmental_vars", file = file.path(data_dir, "germplasm_origin_bioclim_data.RData"))

# Save geographic coordinates with sample names
# Load the germplasm metadata (from the 01_prep_data.R script)
load(file.path(data_dir, "germplasm_metadata.RData"))
load(file.path(data_dir, "gc_marker_data.RData"))


# Create the location file
germ_meta %>%
  left_join(., select(eaa_environmental_data, origin_name, location_abbr)) %>%
  filter(individual %in% row.names(geno_mat_wild)) %>%
  select(individual, origin_latitude, origin_longitude) %>%
  mutate(familyID = 0, paternalID = 0, maternalID = 0, sex = 0, phenotype = -9) %>%
  select(familyID, individualID = individual, paternalID, maternalID, sex, phenotype,
         lat = origin_latitude, long = origin_longitude) %>%
  write_tsv(x = ., file = file.path(data_dir, "wild_cranberry_spa_location_input.txt"), col_names = FALSE)



# Save raster files for each variable -------------------------------------

# Range in lat/long for locations
lat_range <- range(pretty(eaa_environmental_data$latitude))
long_range <- range(pretty(eaa_environmental_data$longitude))

# Build a bounding box for cropping
location_points <- cbind(rep(long_range, each = 2), c(lat_range, rev(lat_range))) %>%
  rbind(., head(., 1)) %>%
  list() %>%
  st_polygon() %>%
  st_sfc(
    crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
  )

# Create a bounding box
location_bbox <-  st_bbox(location_points)

# Convert to igh projection
location_igh <- st_transform(x = location_points, crs = '+proj=igh +lat_0=0 +lon_0=0 +datum=WGS84 +units=m +no_defs')
location_igh_bbox <- st_bbox(location_igh)


# Elevation data
# Read in the elevation data
usa_alt_data <- getData(name = "alt", country = "USA", download = FALSE, keepzip = TRUE,
                        path = file.path(env_dir, "Elevation"))
can_alt_data <- getData(name = "alt", country = "CAN", path = file.path(env_dir, "Elevation"),
                        keepzip = TRUE)


# Combine; usa[[1]] is continental
alt_data <- raster::merge(x = usa_alt_data[[1]], y = can_alt_data)
# Reproject
alt_data1 <- projectRaster(alt_data, crs = '+proj=igh +lat_0=0 +lon_0=0 +datum=WGS84 +units=m +no_defs')

# Crop this
elevation_data_cropped <- crop(x = alt_data1, y = location_igh_bbox)

# WorldClim data
# Read in already downloaded worldclim data
worldclim_dat <- getData(name = "worldclim", download = FALSE, var = "bio", res = "5",
                         path = file.path(env_dir, "WorldClim"))
# Reproject
worldclim_dat1 <- projectRaster(worldclim_dat, crs = '+proj=igh +lat_0=0 +lon_0=0 +datum=WGS84 +units=m +no_defs')

# Crop
worldclim_dat_cropped <- crop(worldclim_dat1, location_igh_bbox)


# Create an empty list to store raster layers
soil_rasters_list <- list()

## Iterate over soil variables
for (vari in soil_vars_df$variable) {

  # Find and read in the raster files for that variable
  raster_files <- list.files(path = file.path(data_dir, "soilgrids"), pattern = vari, full.names = TRUE)
  # get the depth information
  raster_vari_depth <- basename(raster_files) %>%
    str_remove(string = ., pattern = paste0(vari, "_")) %>%
    str_remove(., "_mean.tif")

  raster_list <- lapply(raster_files, raster)
  names(raster_list) <- raster_vari_depth

  # Calculate the interval ranges for weighted means
  vari_depth_intervals <- raster_vari_depth %>%
    str_remove(., "cm") %>%
    map_dbl(~parse(text = .x) %>% eval() %>% abs()) %>%
    setNames(., raster_vari_depth)

  # Aggregate the two raster depths belonging to topsoil and subsoil
  rasters_topsoil <- subset(raster_list, names(raster_list) %in% c("0-5cm", "5-15cm"))
  rasters_subsoil <- subset(raster_list, ! names(raster_list) %in% c("0-5cm", "5-15cm"))

  # Calculate weights
  topsoil_weights <- vari_depth_intervals[names(rasters_topsoil)] %>% {. / sum(.)}
  subsoil_weights <- vari_depth_intervals[names(rasters_subsoil)] %>% {. / sum(.)}

  # Take the weighted average
  rasters_topsoil_avg <- map2(.x = rasters_topsoil, .y = topsoil_weights, `*`) %>%
    reduce(.x = ., .f = `+`)
  rasters_subsoil_avg <- map2(.x = rasters_subsoil, .y = subsoil_weights, `*`) %>%
    reduce(.x = ., .f = `+`)

  # Rename;
  # Add to the list
  names(rasters_topsoil_avg) <- paste0(vari, "_topsoil")
  names(rasters_subsoil_avg) <- paste0(vari, "_subsoil")

  soil_rasters_list[[vari]] <- list(rasters_topsoil_avg, rasters_subsoil_avg)

}


# Merge all soil rasters in to a single raster brick
soilgrids_raster_brick <- unlist(soil_rasters_list) %>%
  `names<-`(., sapply(., names)) %>%
  brick(x = .)
# Rename

raster_dir <- file.path(data_dir, "rasterFiles")
dir.create(path = raster_dir)

# Write everything to disk
for (i in seq_len(nlayers(soilgrids_raster_brick))) {
  writeRaster(x = soilgrids_raster_brick[[i]],
              filename = file.path(raster_dir, paste0("soilgrids_", names(soilgrids_raster_brick)[i], "_data_raster.tif")))

}
for (i in seq_len(nlayers(worldclim_dat_cropped))) {
  writeRaster(x = worldclim_dat_cropped[[i]],
              filename = file.path(raster_dir, paste0("worldclim_", names(worldclim_dat_cropped)[i], "_data_raster.tif")))

}
writeRaster(x = elevation_data_cropped, filename = file.path(raster_dir, "elevation_data_raster.tif"))

