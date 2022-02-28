
<!-- README.md is generated from README.Rmd. Please edit that file -->

# CranberryGermplasmEAA - Scripts

<!-- badges: start -->
<!-- badges: end -->

Below is a description of scripts used in this analysis. Scripts are
listed in order of intended execution.

1.  `01_prep_data.R` - Reads in the starting data from the `Data`
    directory. Filters SNP marker data and collects the population
    metadata.  
2.  `02_get_environmental_data.R` - Queries the WorldClim and SoilGrids
    databases for bioclimatic and soil data for use in the environmental
    association anaylsis.  
3.  `03_population_genetics_statistics.R` - Runs population genetic
    analyses, including *F*<sub>*S**T*</sub> and SPA. Also calculates
    pairwise genetic distance between and within populations.  
4.  `04_lmm_environmental_association.R` - Runs the environmental
    association analysis using the mixed model framework.  
5.  `05_analyze_environmental_association.R` - Performs some analysis of
    the envrionmental association results, including a search for nearby
    annotate genes.

## Figures

`figures.R` - Produces all of the figures and tables (main and
supplemental) that appear in the manuscript.
