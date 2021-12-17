# Germplasm collection environmental association
#
# Conduct the environmental association
#
# Use the R.Sambada package to conduct the environmental association analysis
#


# Load packages
library(raster)
library(tidyverse)
library(readxl)
library(R.SamBada)


# Set directories
source("startup.R")

# Path to sambada binaries
sambada_bin_dir <- file.path(getwd(), "Environment/binaries")

# Directory in which to store sambada inputs
sambada_data_dir <- file.path(data_dir, "sambada_data")


# Load and edit bioclim data as in the LMM script -------------------------

# Read in the marker data
load(file.path(data_dir, "gc_marker_data.RData"))
# Load the germplasm metadata (from the 01_prep_data.R script)
load(file.path(data_dir, "germplasm_metadata.RData"))
# Load the worldclim data
load(file.path(data_dir, "germplasm_origin_bioclim_data.RData"))

# Center and scale
eaa_environmental_data1 <- eaa_environmental_data %>%
  select(location_abbr, all_of(eaa_environmental_vars$variable)) %>%
  mutate(elevation = ifelse(is.na(elevation), 0, elevation))

# Read in the germplasm metadata
germ_meta <- germ_meta %>%
  left_join(., select(eaa_environmental_data, origin_name, location_abbr)) %>%
  # Get the state/province name from the location
  mutate(state = str_sub(string = location_abbr, start = -2, end = -1))

# Prepare the bioclim data as if it were trait data
germplasm_bioclim_data <- germ_meta %>%
  filter(category == "wild", !is.na(origin_latitude)) %>%
  select(individual, latitude = origin_latitude, longitude = origin_longitude, state) %>%
  left_join(., eaa_environmental_data1) %>%
  as.data.frame() %>%
  filter(individual %in% colnames(K_wild))



# Prepare the genomic file ------------------------------------------------

# First change working directories to the data_dir
setwd(dir = data_dir)

# Use the vcf file in the data_dir
geno_file <- list.files(".", "wild_cranberry_cleaned_genotypes.vcf", full.names = TRUE)
# Unzip if gz
if (endsWith(geno_file, ".gz")) {
  R.utils::gunzip(geno_file)
  geno_file <- str_remove(geno_file, ".gz$")
}
# Output file
output_file <- "wild_cranberry_sambada_geno_file.csv"


prepareGeno(fileName = geno_file, outputFile = output_file, saveGDS = TRUE, interactiveChecks = TRUE,
            verbose = TRUE, directory = sambada_bin_dir)



# Move the outputs files --------------------------------------------------

# Move the output to the sambada data dir and change the working directory to that
# directory
file.rename(from = output_file, to = file.path(sambada_data_dir, output_file))
gds_file <- list.files(path = ".", pattern = ".gds$", full.names = TRUE)
file.rename(from = gds_file, to = file.path(sambada_data_dir, gds_file))

setwd(sambada_data_dir)


# Prepare the environmental data ------------------------------------------

# This code will essentially get around the createEnv function for assembling a
# file with the environmental data. Here we will create a space-delimited file with
# the location names and environmental data
env_file <- "./wild_cranberry_sambada_env_file.csv"
gds_file <- list.files(".", ".gds$", full.names = TRUE)

germplasm_bioclim_data %>%
  select(ID_indiv = individual, latitude, longitude, all_of(eaa_environmental_vars$variable)) %>%
  distinct() %>%
  write_delim(x = ., file = env_file, delim = " ")

# Output file
output_file <- "wild_cranberry_sambada_env_file_export.csv"

# Now run prepareEnv
prepareEnv(envFile = env_file, outputFile = output_file, maxCorr = 0.85, idName = "ID_indiv",
           genoFile = gds_file, numPc = 0.2, x = "latitude", y = "longitude", locationProj=4326,
           mafThresh = 2 / n_distinct(germplasm_bioclim_data$individual), numPop = NULL, missingnessThresh = 0.1,
           interactiveChecks = TRUE, verbose = TRUE)




# Run sambada -------------------------------------------------------------

# Reset the geno file
geno_file <- "wild_cranberry_sambada_geno_file.csv"
env_file <- output_file
output_file <- "wild_cranberry_sambada_results"

sambadaParallel(genoFile = geno_file, envFile = env_file, idGeno = "ID_indiv", idEnv = "ID_indiv",
                outputFile = output_file, dimMax = 2, saveType = "END ALL", populationVar = "LAST",
                directory = sambada_bin_dir)




# Prepare the output ------------------------------------------------------

sambada_output <- prepareOutput(sambadaname = output_file, dimMax = 2, gdsFile = gds_file,
                                popStr = TRUE, interactiveChecks = FALSE)

# Get the data output
sambada_output_markers <- sambada_output$sambadaOutput %>%
  as_tibble() %>%
  rename(variable = Env_1, marker = snp, marker_allele = Marker)

# Variables with output
sort(unique(sambada_output_markers$variable))

sambada_output_markers %>%
  filter(variable == "bio15") %>%
  rename(Chrom = chr, Position = pos) %>%
  mutate(p.val = -log10(qvalueG)) %>%
  sommer::manhattan(map = .)
















