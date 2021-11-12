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


# Get biophysical soil variables ------------------------------------------


# Load additional libraries
library(XML)
library(rgdal)
library(gdalUtils)
library(sf)

# Set the gdal installation location
gdal_setInstallation(search_path = "C:\\Users\\jeffrey.neyhart\\Anaconda3\\Library\\bin", verbose = TRUE)

# URL for soil variables
url <- "https://files.isric.org/soilgrids/latest/data/"
# Projection
igh <- '+proj=igh +lat_0=0 +lon_0=0 +datum=WGS84 +units=m +no_defs'




# List of subdirectories of interest
soil_vars <- c(
  "bdod" = "bulk_density",
  "cec" = "cation_exchange_cap",
  "cfvo" = "course_fragment_content",
  "clay" = "clay",
  "nitrogen" = "total_nitrogen",
  "ocd" = "organic_carbon_density",
  "phh2o" = "pH",
  "sand" = "sand",
  "silt" = "silt",
  "soc" = "soil_organic_carbon"
)

soil_vars_df <- tibble(variable = names(soil_vars), full_name = soil_vars, class = "soil")

# Get data from each layer
depths <- c("0-5cm", "5-15cm", "15-30cm", "30-60cm", "60-100cm", "100-200cm")

# quantiles of interest
quantiles <- c("mean")

# Combinations
soil_var_combs <- crossing(soil_vars_df, depth = depths, quantile = quantiles)

## Download files - only do this once ##

# Create a new storage directory
soil_data_dir <- file.path(data_dir, "soilgrids")
dir.create(path = soil_data_dir)

# Iterate over combinations
for (i in seq_len(nrow(soil_var_combs))) {
  # Get the variable of interest
  voi = soil_var_combs$variable[i] # variable of interest
  depth = soil_var_combs$depth[i]
  quantile = soil_var_combs$quantile[i]

  voi_layer = paste(voi,depth,quantile, sep="_")
  print(voi_layer)

  wcs_path = paste0("https://maps.isric.org/mapserv?map=/map/",voi,".map") # Path to the WCS. See maps.isric.org
  wcs_service = "SERVICE=WCS"
  wcs_version = "VERSION=2.0.1"
  wcs = paste(wcs_path,wcs_service,wcs_version,sep="&") # This works for gdal >= 2.3

  l1 <- newXMLNode("WCS_GDAL")
  l1.s <- newXMLNode("ServiceURL", wcs, parent=l1)
  l1.l <- newXMLNode("CoverageName", voi_layer, parent=l1)


  # Save to local disk
  xml.out = file.path(soil_data_dir, "sg.xml")
  saveXML(l1, file = xml.out)


  long_range <- range(pretty(location_metadata$longitude))
  lat_range <- range(pretty(location_metadata$latitude))

  ## Build the bounding box
  location_points <- cbind(rep(long_range, each = 2), c(lat_range, rev(lat_range))) %>%
    rbind(., head(., 1)) %>%
    list() %>%
    st_polygon() %>%
    st_sfc(
      crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
    )

  location_igh <- st_transform(x = location_points, crs = '+proj=igh +lat_0=0 +lon_0=0 +datum=WGS84 +units=m +no_defs')

  location_bbox <- st_bbox(obj = location_igh)
  # Create the query bbox
  ulx = location_bbox$xmin
  uly = location_bbox$ymax
  lrx= location_bbox$xmax
  lry = location_bbox$ymin
  bb <- c(ulx, uly, lrx, lry)

  # Get the whole region as a geotiff
  # Download raster as GeoTIFF (Warning: it can be large!)
  file.out <- file.path(soil_data_dir, paste0(voi_layer, ".tif"))
  if (file.exists(file.out)) next

  gdal_translate(xml.out, file.out,
                 tr=c(250,250), projwin=bb,
                 projwin_srs = igh, co=c("TILED=YES","COMPRESS=DEFLATE","PREDICTOR=2","BIGTIFF=YES"),
                 verbose=TRUE )

}

## Extract data from the raster files
# List the raster files
raster_files <- list.files(path = soil_data_dir, full.names = TRUE, pattern = ".tif")

# Iterate over locations and subset the raster
# Subset data for each coordinate
location_metadata1 <- location_metadata %>%
  filter_at(vars(latitude, longitude), all_vars(!is.na(.))) %>%
  mutate(soil_data = list(NULL))

pb <- progress::progress_bar$new(total = nrow(location_metadata1))
# Iterate over rows
for (i in seq_len(nrow(location_metadata1))) {

  # Location code
  loc_code_i <- location_metadata1$location_abbr[i]

  # Get latitude/longitude
  lat_i <- location_metadata1$latitude[i]
  long_i <- location_metadata1$longitude[i]

  # Build extent (0.1 degrees)
  geo_add <- (soil_data_aggr_res / 60) / 2
  long_range <- c(long_i - geo_add, long_i + geo_add)
  lat_range <- c(lat_i - geo_add, lat_i + geo_add)

  ## Build the bounding box
  location_points <- cbind(rep(long_range, each = 2), c(lat_range, rev(lat_range))) %>%
    rbind(., head(., 1)) %>%
    list() %>%
    st_polygon() %>%
    st_sfc(
      crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
    )

  location_igh <- st_transform(x = location_points, crs = '+proj=igh +lat_0=0 +lon_0=0 +datum=WGS84 +units=m +no_defs')
  location_bbox <- st_bbox(obj = location_igh)


  # Output df of raster data
  raster_data_df_i <- str_remove(basename(raster_files), ".tif") %>%
    str_split(string = ., pattern = "_") %>%
    do.call("rbind", .) %>%
    as.data.frame() %>%
    `names<-`(., c("variable", "depth", "quantile")) %>%
    mutate(value = as.numeric(NA))

  # Iterate over the raster files
  for (j in seq_along(raster_files)) {

    # Read in the raster file
    raster_in <- raster(raster_files[j])

    # Use the coordinates of the location to subset the raster
    raster_subset <- crop(x = raster_in, y = location_bbox)
    # Take the mean of non-zero values
    raster_mean_value <- values(raster_subset) %>%
      subset(. > 0) %>%
      mean

    # add this to the df
    raster_data_df_i$value[j] <- raster_mean_value

  }

  # Add the df to the larger df
  location_metadata1$soil_data[[i]] <- raster_data_df_i

  pb$tick()

}

# Unnest the soil data
location_soil_data <- unnest(location_metadata1) %>%
  # Calculate the range of each interval
  mutate(depth_range = str_remove(depth, "cm"),
         depth_range = map_dbl(depth_range, ~abs(eval(parse(text = .x))))) %>%
  spread(variable, value) %>%
  # Convert to standard units, then
  mutate_at(vars(bdod, nitrogen), ~ . / 100) %>%
  mutate_at(vars(cec, cfvo, clay, phh2o, sand, silt, soc, ocd), ~. / 10)

# aggregate by interval in topsoil and subsoil
location_soil_data_aggr <- location_soil_data %>%
  mutate(soil_layer = ifelse(depth %in% c("0-5cm", "5-15cm"), "topsoil", "subsoil")) %>%
  group_by(location_abbr, soil_layer) %>%
  summarize_at(vars(all_of(names(soil_vars))), ~weighted.mean(x = ., w = depth_range, na.rm = TRUE)) %>%
  ungroup()




# Save this
save("location_soil_data", "location_soil_data_aggr", "soil_vars_df",
     file = file.path(data_dir, "germplasm_origin_soil_data.RData"))








