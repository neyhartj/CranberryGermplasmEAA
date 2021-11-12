#!/bin/bash

## Wild cranberry SPatial Ancestry Analysis (SPA)
##
## This script runs the SPA software using the genotypes of the wild cranberries
## along with their known geographic origins
##
##

# Path to the input files
pedbase="Data/wild_cranberry_cleaned_genotypes"

# Path to the location file
locfile="Data/wild_cranberry_spa_location_input.txt"

# Path to the SPA software
spa="../../../../../Software/spa/spa.exe"


# Run spa
$spa --pfile $pedbase --location-input $locfile --model-output Results/wild_cranberry_knownGeo_spa.out -v 1

# Run spa without known geographic origins
$spa --pfile $pedbase --location-output Results/wild_cranberry_inferGeno_locInfo_spa.out --model-output Results/wild_cranberry_inferGeo_spa.out -v 1 -r 1e-6
