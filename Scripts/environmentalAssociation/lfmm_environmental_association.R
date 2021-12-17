# Germplasm collection environmental association
#
# Conduct the environmental association
#
# Use the linear mixed model and latent factor mixed model frameworks for
# environmental association analysis
#


# Load packages
library(tidyverse)
library(readxl)
library(vcfR)
library(LEA)
library(snps)
library(R.SamBada)
library(qvalue)

# Load the startup script
source("startup.R")



# Read in the marker data
load(file.path(data_dir, "gc_marker_data.RData"))

# Load the worldclim data
load(file.path(data_dir, "germplasm_origin_bioclim_data.RData"))

# Center and scale
eaa_environmental_data1 <- eaa_environmental_data %>%
  select(location_abbr, elevation, latitude, longitude, contains("bio"), starts_with("IC")) %>%
  mutate_at(vars(contains("bio")), scale, scale = FALSE) %>%
  mutate_at(vars(contains("bio")), as.numeric) %>%
  mutate(elevation = ifelse(is.na(elevation), 0, elevation))

# Read in the germplasm metadata
germ_meta <- read_csv(file = file.path(data_dir, "gc_metadata.csv")) %>%
  left_join(., select(eaa_environmental_data, origin_name, location_abbr))

# Read in variety origin information
native_sel_meta <- read_csv(file = file.path(data_dir, "native_selection_metadata.csv"))


# Prepare the bioclim data ------------------------------------------------

# Prepare the bioclim data as if it were trait data
germplasm_bioclim_data <- germ_meta %>%
  filter(category == "wild", !is.na(origin_latitude)) %>%
  select(individual, latitude = origin_latitude, longitude = origin_longitude) %>%
  left_join(., eaa_environmental_data1) %>%
  as.data.frame() %>%
  filter(individual %in% colnames(K))

# Manipulate the marker data -------------------------------------------------

geno_mat_wild1 <- geno_mat_wild[germplasm_bioclim_data$individual,]

geno_hmp_wild1 <- geno_hmp_wild %>%
  select(marker:alleles, row.names(geno_mat_wild1))

# Split the marker matrix by chromosome
geno_mat_list <- geno_hmp_wild %>%
  select(marker, chrom) %>%
  split(.$chrom) %>% map("marker") %>%
  map(~geno_mat_wild1[,.x, drop = FALSE])

# Subset the geno_hmp_wild object
geno <- geno_hmp_wild %>%
  select(marker, chrom, pos, all_of(colnames(K))) %>%
  mutate_at(vars(colnames(K)), ~. - 1) %>%
  as.data.frame()

# Subset the native selections
geno_mat_native <- geno_mat_all[native_sel_meta$germplasm_collection_name,]






# Run the latent factor mixed model analysis ------------------------------

## Format input data ##

# Convert the marker matrix to 0, 1, 2
geno_mat_wild2 <- geno_mat_wild1 + 1

# Save as a space-separated txt file
lfmm_file <- write.lfmm(R = geno_mat_wild2, output.file = file.path(data_dir, "lfmm_marker_genotypes.lfmm"))

# Convert the lfmm data to geno data
lfmm_genotype <- lfmm2geno(input.file = file.path(data_dir, "lfmm_marker_genotypes.lfmm"))

# Convert the bioclim variables in to a matrix
bioclim_X <- germplasm_bioclim_data %>%
  select(individual, all_of(eaa_environmental_vars$variable)) %>%
  column_to_rownames("individual") %>%
  as.matrix()

# Save as .env file
write.env(R = bioclim_X[,!startsWith(colnames(bioclim_X), "IC")], output.file = file.path(data_dir, "lfmm_eco_data.env"))
# Just the linear combinations
write.env(R = bioclim_X[,startsWith(colnames(bioclim_X), "IC")], output.file = file.path(data_dir, "lfmm_eco_data_IC.env"))


## Analyze population structure
pc <- pca(input.file = lfmm_file, scale = TRUE)

# Use tracy-widom test for eigenvalues
## Results suggest K = 13

# Create a class snmf
obj_snmf = snmf(input.file = lfmm_genotype, K = seq_len(n_distinct(germplasm_bioclim_data$location_abbr)),
                entropy = TRUE, ploidy = 2, project="new")

# Plot the results
plot(obj_snmf)
# It appears that K = 7 is appropriate
K_min_entropy <- which.min(map_dbl(obj_snmf@runs, slot, "crossEntropy"))
# Visualize ancestry
barplot(t(Q(obj_snmf, K = K_min_entropy)), col = 1:6)

## Population differentiation test
p_diff <- snmf.pvalues(object = obj_snmf, entropy = TRUE, ploidy = 2, K = K_min_entropy)
pvalues = p_diff$pvalues
plot(-log10(pvalues), pch = 19, col = "blue", cex = .5)




# Run the LFMM function




### Alternative using the lfmm2 function ###


# Convert the bioclim variables in to a matrix
bioclim_X <- germplasm_bioclim_data %>%
  select(individual, all_of(biovars)) %>%
  column_to_rownames("individual") %>%
  as.matrix()

bioclim_X1 <- bioclim_X[,-c(17:19),drop = FALSE]

# Iterate over environmental variables
#
# First create a list
lfmm_out_list <- list()

for (i in seq_len(ncol(bioclim_X1))) {
  # Extract the environmental variable
  X <- bioclim_X1[,i, drop = FALSE]

  # Fit the latent factor mixed model
  lfmm_out <- lfmm2(input = geno_mat_wild2, env = X, K = K_min_entropy)

  # Comput p-values
  lfmm_test_out <- lfmm2.test(object = lfmm_out, input = geno_mat_wild2, env = X,
                              genomic.control = TRUE, linear = TRUE)

  # Add the results to the list
  lfmm_out_list[[i]] <- lfmm_test_out

}

# Add names
names(lfmm_out_list) <- colnames(bioclim_X1)

# Get the genomic inflation factor
map_dbl(lfmm_out_list, "gif")


# Plot
for (vari in names(lfmm_out_list)) {

  # p-values
  pvals <- lfmm_out_list[[vari]]$pvalues

  # combine pvalues with the map
  scores1 <- as.data.frame(pvals) %>%
    rownames_to_column("marker") %>%
    rename(pvalue = "V1") %>%
    mutate(marker = str_remove(marker, "Response "),
           score = -log10(pvalue)) %>%
    left_join(., select(geno_hmp_wild1, marker, Chrom = chrom, Position = pos), by = "marker")

  # Manahattan plot using sommer
  manhattan(map = scores1, PVCN = "score", main = bioclim_vars[vari], cex = 0.5)

}
















# Run the r.sambada analysis ----------------------------------------------

# Create a temporary directory to store files
temp_dir <- "C:/Users/jeffrey.neyhart/Documents/TemporaryData/"


## Prep Data ##

## Prep the genomic data ##
# Read in the VCF file
vcf_in <- read.vcfR(file = file.path(geno_dir, "phased_marker_genotype_db.vcf"))

# Subset the markers and individuals from the environmental GWAS
rows <- which(vcf_in@fix[,"ID"] %in% colnames(geno_mat_wild1))
cols <- c(1, which(colnames(vcf_in@gt) %in% row.names(geno_mat_wild1)))

vcf_in1 <- vcf_in[rows, cols]

vcf_in1@fix[,"CHROM"] <- as.character(as.numeric(vcf_in1@fix[,"CHROM"]))

# Save this VCF
write.vcf(x = vcf_in1, file = file.path(temp_dir, "gcWild_filtered_phased_markers.vcf"))


# Use this vcf to create a csv file for use in sambada
# This file will be space-delimited with the first column being the individual name
# and subsequent columns designated SNPs as triplets: one column for each of the 3 homo/het
# genotype states; elements of the file are 0 (for false) or 1 (for true) designating
# whether the jth individual carries that genotype at the ith SNP

# Tranpose the genotype matrix
gtT <- t(vcf_in1@gt[,-1])
# Get the SNP metdata
snp_meta <- vcf_in1@fix[,c("ID", "REF", "ALT")]

# Create a empty matrix to store the results
geno_output_list <- replicate(n = nrow(snp_meta), matrix(data = 0, nrow = nrow(gtT), ncol = 3),
                              simplify = FALSE)
geno_output_snp_allele_names_list <- list()

# Iterate over SNPs
for (i in seq_len(nrow(snp_meta))) {
  # Get the SNP name
  snp_name <- snp_meta[i,"ID"]
  # Create names with the ref/alt combinations
  alleles <- snp_meta[i,c("REF", "ALT")]
  allele_suf <- c(paste0(alleles[c(2,2)], collapse = ""), paste0(alleles[c(1,2)], collapse = ""),
                  paste0(alleles[c(1,1)], collapse = ""))
  snp_allele_names <- paste0(snp_name, "_", allele_suf)

  # Create the incidence matrix
  snp_calls <- sapply(X = strsplit(x = gtT[,i], split = "\\|"), FUN = function(snp) sum(as.numeric(snp))) + 1
  for (j in seq_along(snp_calls)) {
    geno_output_list[[i]][j,snp_calls[j]] <- 1
  }

  # add the colnames to a list
  geno_output_snp_allele_names_list[[i]] <- snp_allele_names

}

# cbind
geno_output_mat1 <- do.call("cbind", geno_output_list)

# Edit the dimnames
dimnames(geno_output_mat1) <- list(row.names(gtT), unlist(geno_output_snp_allele_names_list))

# Save this to a csv
geno_output_mat1 %>%
  as.data.frame() %>%
  rownames_to_column("individual") %>%
  select(1:1000) %>%
  write_delim(x = ., path = file.path(temp_dir, "gcWild_genos_prepared_manual.csv"), delim = " ")




# Read the genomic data into sambada
prepareGeno(fileName = file.path(temp_dir, "gcWild_filtered_phased_markers.vcf"),
            outputFile = file.path(temp_dir, "gcWild_genos_prepared.csv"), saveGDS = TRUE,
            mafThresh = 0, missingnessThresh = 1, interactiveChecks = FALSE,
            directory = file.path(proj_dir, "Environment/binaries"), verbose = FALSE)

# This is the true output file name (since it saved the GDS file)
gds_output <- file.path(temp_dir, "gcWild_filtered_phased_markers.gds")

# Move the file
file.copy(from = file.path(getwd(), "gcWild_filtered_phased_markers.gds"),
          to = gds_output, overwrite = TRUE)
# Remove the original
unlink(x = file.path(getwd(), "gcWild_filtered_phased_markers.gds"))

## Copy other files
file.copy(from = file.path(tempdir(), "gcWild_filtered_phased_markers_filtered.ped"),
          to = file.path(temp_dir, "gcWild_filtered_phased_markers_filtered.ped"), overwrite = T)
file.copy(from = file.path(tempdir(), "gcWild_filtered_phased_markers_filtered.map"),
          to = file.path(temp_dir, "gcWild_filtered_phased_markers_filtered.map"), overwrite = T)

# Rerun prepareGeno with the ped file
prepareGeno(fileName = file.path(temp_dir, "gcWild_filtered_phased_markers_filtered.ped"),
            outputFile = file.path(temp_dir, "gcWild_genos_prepared.csv"), saveGDS = FALSE,
            mafThresh = 0, missingnessThresh = 1, interactiveChecks = FALSE,
            directory = file.path(proj_dir, "Environment/binaries"), verbose = FALSE)



## Prep the environmental data ##

# Output an environment file with the bioclim data; this will be space-delimited
# where rows are individuals and columns are bioclim variable values for the originating
# location of those individuals
germplasm_bioclim_data1 %>%
  select(individual, latitude, longitude, contains("bio")) %>%
  write_delim(x = ., path = file.path(temp_dir, "location_file_env.csv"), delim = " ")


# Run the preparEnv function
prepareEnv(envFile = file.path(temp_dir, "location_file_env.csv"),
           outputFile = file.path(temp_dir, "location_file_env_export.csv"),
           genoFile = gds_output, idName = "individual", separator = " ",
           maxCorr = 0.85, numPc = 0.2, numPop=NULL,
           mafThresh = 0, missingnessThresh = 1, ldThresh = 1,
           x='longitude', y='latitude',
           interactiveChecks=FALSE, locationProj=4326 )

# Run sambada
sambadaParallel(genoFile = file.path(temp_dir, "gcWild_genos_prepared.csv"),
                envFile = file.path(temp_dir, "location_file_env_export.csv"),
                outputFile = file.path(temp_dir, "gc_sambada"),
                idGeno='ID_indiv', idEnv='individual', dimMax = 2, cores = 1,
                saveType = "END ALL", populationVar = "LAST", keepAllFiles = TRUE,
                directory = file.path(proj_dir, "Environment/binaries"))


## Copy the output to results
result_files <- list.files(temp_dir, pattern = "gcWild_genos_prepared-", full.names = TRUE)
for (file in result_files) file.copy(from = file, to = file.path(result_dir, basename(file)), overwrite = T)

# Prepare the output
prep <- prepareOutput(sambadaname = file.path(result_dir, "gcWild_genos_prepared"),
                      gdsFile = gds_output, dimMax = 2, popStr = TRUE, interactiveChecks = FALSE)

# Copy
prep_filtered <- prep
prep_filtered$sambadaOutput <- prep_filtered$sambadaOutput %>%
  filter(Gscore > 6)

# Manhattan plots
plotManhattan(preparedOutput = prep, varEnv = sort(unique(prep$sambadaOutput$Env_1)),
              chromo='all',valueName='pvalueW', threshold = 0.05)













