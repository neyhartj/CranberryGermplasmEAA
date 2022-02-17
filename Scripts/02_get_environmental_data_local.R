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

env_dir <- file.path(cran_dir, "EnvironmentalData")

# Resolution at which to aggregate soil data (in minutes)
soil_data_aggr_res <- 2.5

# Exclude monthly variables?
exclude_monthly_bioclim <- TRUE



# Read in data ------------------------------------------------------------

# Read in population metadata - local
load(file.path(data_dir, "population_metadata_and_genotypes.RData"))


# Create abbreviations for the origin names
location_metadata <- pop_metadata %>%
  # filter for missing origin or coordinates; select only genotypes in K
  filter(!is.na(location_of_origin),
         !is.na(latitude),
         individual %in% row.names(K_wild)) %>%
  distinct(location_of_origin, latitude, longitude) %>%
  mutate(location_of_origin1 = str_remove(location_of_origin, "west branch, ")) %>%
  separate(location_of_origin1, c("location", "state"), sep = ", ") %>%
  mutate(origin_abbr = abbreviate(location, 3),
         origin_abbr = paste0(toupper(origin_abbr), state)) %>%
  select(location, location_of_origin, state, location_abbr = origin_abbr, latitude, longitude)

# Set the min max of lat/long for subsetting data
lat_range <- range(location_metadata$latitude, na.rm = TRUE)
long_range <- range(location_metadata$longitude, na.rm = TRUE)



# Format environmental variables ------------------------------------------

# Names of bioclim variables
environ_vars_tbl <- tribble(
  ~ variable, ~ full_name, ~ unit,
  "elevation", "Elevation", "m",
  "latitude", "Latitude", "",
  "longitude", "Longitude", "",
  "bio1", 'Annual Mean Temperature', "C",
  "bio2", 'Mean Diurnal Range', "C",
  "bio3", 'Isothermality', "",
  "bio4", 'Temperature Seasonality', "C x 100",
  "bio5", 'Max Temperature of Warmest Month',  "C",
  "bio6", 'Min Temperature of Coldest Month', "C",
  "bio7", 'Temperature Annual Range', "C",
  "bio8", 'Mean Temperature of Wettest Quarter', "C",
  "bio9", 'Mean Temperature of Driest Quarter', "C",
  "bio10", 'Mean Temperature of Warmest Quarter', "C",
  "bio11", 'Mean Temperature of Coldest Quarter', "C",
  "bio12", 'Annual Precipitation', "mm",
  "bio13", 'Precipitation of Wettest Month', "mm",
  "bio14", 'Precipitation of Driest Month', "mm",
  "bio15", 'Precipitation Seasonality', "mm",
  "bio16", 'Precipitation of Wettest Quarter', "mm",
  "bio17", 'Precipitation of Driest Quarter', "mm",
  "bio18", 'Precipitation of Warmest Quarter', "mm",
  "bio19", 'Precipitation of Coldest Quarter', "mm"
)

environ_vars <- setNames(environ_vars_tbl$full_name, environ_vars_tbl$variable)


# Soil variables
soil_vars_tbl <- tribble(
  ~ variable, ~ full_name, ~ unit,
  "bdod", "bulk_density", "cg cm^-3",
  "cec", "cation_exchange_cap", "mmol(c) kg^-1",
  "cfvo", "course_fragment_content", "cm^3 dm^-3",
  "clay", "clay", "g kg^-1",
  "nitrogen", "total_nitrogen", "cg kg^-1",
  "ocd", "organic_carbon_density", "hg m^-3",
  "phh2o", "pH", "pH x 10",
  "sand", "sand", "g kg^-1",
  "silt", "silt", "g kg^-1",
  "soc", "soil_organic_carbon", "dg kg^-1"
)


soil_vars <- setNames(soil_vars_tbl$full_name, soil_vars_tbl$variable)




# Get altitude and WorldClim data ------------------------------------------------------

# Many lines of this code are from here:
# https://github.com/MorrellLAB/Env_Assoc/blob/master/script/GWAS/worldclim.miner3.r

# Get worldclim data, bioclimatic variables at 5 minute resolution (~9 sq km)
#
#
# # This is only to be used once; Use the next command to read in data
# worldclim_dat <- getData(name = "worldclim", var = "bio", res = "5",
#                          path = file.path(env_dir, "WorldClim"))
#
# # Tmin, tmax, precip
# worldclim_dat <- getData(name = "worldclim", var = "tmin", res = "5",
#                          path = file.path(env_dir, "WorldClim"))
#
# worldclim_dat <- getData(name = "worldclim", var = "tmean", res = "5",
#                          path = file.path(env_dir, "WorldClim"))
#
# worldclim_dat <- getData(name = "worldclim", var = "tmax", res = "5",
#                          path = file.path(env_dir, "WorldClim"))
#
# worldclim_dat <- getData(name = "worldclim", var = "prec", res = "5",
#                          path = file.path(env_dir, "WorldClim"))

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

worldclim_dat_tmin <- getData(name = "worldclim", var = "tmin", res = "5", download = FALSE,
                              path = file.path(env_dir, "WorldClim"))

worldclim_dat_tmean <- getData(name = "worldclim", var = "tmean", res = "5", download = FALSE,
                               path = file.path(env_dir, "WorldClim"))

worldclim_dat_tmax <- getData(name = "worldclim", var = "tmax", res = "5", download = FALSE,
                              path = file.path(env_dir, "WorldClim"))

worldclim_dat_prec <- getData(name = "worldclim", var = "prec", res = "5", download = FALSE,
                              path = file.path(env_dir, "WorldClim"))


# Combine all raster stacks
worldclim_dat1 <- stack(worldclim_dat, worldclim_dat_tmin, worldclim_dat_tmean, worldclim_dat_tmax,
                        worldclim_dat_prec)


# Crop this data
long_range1 <- range(pretty(long_range))
lat_range1 <- range(pretty(lat_range))
extent1 <- extent(list(x = long_range1, y = lat_range1))
worldclim_dat1 <- crop(x = worldclim_dat1, y = extent1)
alt_data1 <- crop(x = alt_data, y = extent1)




# Combine with temp and precip
environ_vars <- c(environ_vars,
                  setNames(paste0(month.abb, " Min Temperature"), names(worldclim_dat_tmin)),
                  setNames(paste0(month.abb, " Max Temperature"), names(worldclim_dat_tmax)),
                  setNames(paste0(month.abb, " Mean Temperature"), names(worldclim_dat_tmean)),
                  setNames(paste0(month.abb, " Precipitation"), names(worldclim_dat_prec))
                  )


# Combine into tidy df with classifiers for each variable
environ_vars_df <- tibble(variable = names(environ_vars), full_name = environ_vars) %>%
  mutate(class = case_when(variable %in% paste0("bio", 1:11) ~ "temperature",
                           variable %in% paste0("bio", 12:19) ~ "precipitation",
                           str_detect(full_name, "Temperature") ~ "temperature",
                           str_detect(full_name, "Precipitation") ~ "precipitation",
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



# Get soilgrids data ------------------------------------------------------

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



soil_vars_df <- tibble(variable = names(soil_vars), full_name = soil_vars, class = "soil")

# Get data from each layer
depths <- c("0-5cm", "5-15cm", "15-30cm", "30-60cm", "60-100cm", "100-200cm")

# quantiles of interest
quantiles <- c("mean")

# Combinations
soil_var_combs <- crossing(soil_vars_df, depth = depths, quantile = quantiles)

## Download files - only do this once ##

# Create a new storage directory
soil_data_dir <- file.path(env_dir, "SoilGrids")
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




# Process and combine the worldclim and soilgrids  ------------------------

library(sf)


# Spread out the aggregated soil data
location_soil_data_aggr1 <- location_soil_data_aggr %>%
  gather(variable, value, soil_vars_df$variable) %>%
  unite(var, variable, soil_layer, sep = "_") %>%
  spread(var, value)

# Combine the worldclim and soil data
location_bioclim_data <- full_join(location_worldclim_data, location_soil_data_aggr1) %>%
  # Elevation is zero if NA
  mutate(elevation = ifelse(is.na(elevation), 0, elevation))

# Edit the soil variable df
soil_vars_df1 <- soil_vars_df %>%
  crossing(., distinct(location_soil_data_aggr, soil_layer)) %>%
  mutate(variable = paste0(variable, "_", soil_layer),
         full_name = paste0(full_name, "_", soil_layer),
         full_name = str_to_title(str_replace_all(full_name, "_", " "))) %>%
  select(-soil_layer)

# Combine the list of variable names
bioclim_vars_df <- bind_rows(environ_vars_df, soil_vars_df1)

# Pull out variables of interest
if (exclude_monthly_bioclim) {
  # Variables of interest for analysis
  env_variables_of_interest <- bioclim_vars_df %>%
    subset(str_detect(full_name, paste0(month.abb, collapse = "|"), negate = TRUE), variable, drop = TRUE)

} else {
  env_variables_of_interest <- bioclim_vars_df$variable

}



# Principal component analysis of all variables ---------------------------

# Center and scale
# Convert the data to a matrix
location_bioclim_data_scaled <- location_bioclim_data %>%
  select(location_abbr, all_of(env_variables_of_interest)) %>%
  # Center and scale the variables
  mutate_at(vars(-location_abbr), ~as.numeric(scale(.)))

location_bioclim_data_scaled_mat <- location_bioclim_data_scaled %>%
  as.data.frame() %>%
  column_to_rownames("location_abbr") %>%
  as.matrix()

## Run PCA
pca_env1 <- prcomp(x = location_bioclim_data_scaled_mat, retx = TRUE)
# Examine proportion of variance explained
pca_env1_propvar <- (pca_env1$sdev^2)/sum(pca_env1$sdev^2)
# Cumulative sum
(pc_varexp <- cumsum(pca_env1_propvar))

# Extract the PCs that together explain more than 90% of the variance
# Tidy
pca_env1_tidy <- broom::tidy(pca_env1) %>%
  filter(PC %in% seq_len(max(which(pc_varexp <= 0.90)) + 1)) %>%
  rename(location_abbr = row)


# Create a table of loadings for each of these PCs
pca_env1_loadings_tidy <- pca_env1$rotation[,unique(pca_env1_tidy$PC),drop = FALSE] %>%
  as.data.frame() %>%
  rownames_to_column("variable") %>%
  gather(PC, loading, -variable) %>%
  as_tibble() %>%
  # arrange by PC and then loading
  arrange(PC, desc(abs(loading)))

# Create a nice table of the variables with the highest-magnitude loading for each PC
pca_env_max_loadings <- pca_env1_loadings_tidy %>%
  group_by(PC) %>%
  top_n(x = ., n = 1, wt = abs(loading)) %>%
  ungroup() %>%
  left_join(., select(eaa_environmental_vars, variable, full_name)) %>%
  mutate(prop_var_exp = pca_env1_propvar[seq_len(nrow(.))]) %>%
  select(PC, prop_var_exp, variable, full_name, loading)


# Append the variables of interest object
env_variables_of_interest <- union(env_variables_of_interest, paste0("PC", unique(pca_env1_tidy$PC)))


# Add the PCs to the location_bioclim_data
eaa_environmental_data <- pca_env1_tidy %>%
  mutate(PC = paste0("PC", PC)) %>%
  spread(PC, value) %>%
  left_join(location_bioclim_data, .) %>%
  select(location:location_abbr, all_of(env_variables_of_interest))

eaa_environmental_vars <- bioclim_vars_df %>%
  add_row(variable = str_subset(colnames(eaa_environmental_data), "PC"),
          full_name = variable, class = "principal_component") %>%
  filter(variable %in% env_variables_of_interest)

# Add units from above
eaa_environmental_vars <- eaa_environmental_vars %>%
  mutate(variable_join = str_remove_all(variable, "_topsoil|_subsoil")) %>%
  left_join(., select(bind_rows(environ_vars_tbl, soil_vars_tbl), variable, unit),
            by = c("variable_join" = "variable")) %>%
  mutate(unit = case_when(
    variable %in% str_subset(names(environ_vars), "^tm") ~ "C",
    variable %in% str_subset(names(environ_vars), "^prec") ~ "mm",
    unit == "" ~ as.character(NA),
    TRUE ~ as.character(unit))
  ) %>%
  select(-variable_join)


# Save everything
save("eaa_environmental_data", "eaa_environmental_vars", "env_variables_of_interest", "pca_env1_loadings_tidy",
     "pca_env_max_loadings",
     file = file.path(data_dir, "germplasm_origin_bioclim_data.RData"))



# Other unused code -------------------------------------





