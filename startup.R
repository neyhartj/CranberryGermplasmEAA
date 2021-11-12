## Project name here
##
## This is a startup script to load relevant data or functions
##


## Packages to load


## Directories
proj_dir <- here::here()


# Other directories
fig_dir <- file.path(proj_dir, "Figures")
data_dir <- file.path(proj_dir, "Data")
result_dir <- file.path(proj_dir, "Results")

# Load functions
source("functions.R")


## Put anything you want below ##


# Cranberry dir
cran_dir <- strsplit(proj_dir, "/")[[1]] %>%
  {.[seq_len(which(. == "CranberryLab"))]} %>%
  paste0(collapse = "/")

# Directory of genotypic data
geno_dir <- file.path(cran_dir, "Genotyping/MarkerDatabase")
