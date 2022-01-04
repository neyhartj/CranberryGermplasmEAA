# Germplasm collection environmental association
#
# Data preparation
#
# 1. Copy phenotypic data
# 2. Copy germplasm collection metadata
#

# Load packages
library(tidyverse)
library(readxl)
library(snps)
library(rrBLUP)
library(vcfR)

# Project directory
proj_dir <- here::here()
data_dir <- file.path(proj_dir, "data")

# CranberryLab directory
path_split <- str_split(string = proj_dir, pattern = "/")[[1]]
cran_dir <- paste0(path_split[seq_len(str_which(string = path_split, pattern = "CranberryLab"))], collapse = "/")

# Directory of genotypic data
geno_dir <- file.path(cran_dir, "Genotyping/MarkerDatabase")


## User parameters
max_LD_r2 <- 0.99999 # Max LD of any one pair of markers
max_LD_r2_2 <- 0.90
min_mac <- 10 # Minimum minor allele count
max_indiv_miss <- 0.7
max_snp_miss <- 0.7


# Load the population metadata --------------------------------------------

pop_metadata <- read_csv(file = file.path(data_dir, "population_metadata.csv"))


# Copy marker genotype data -----------------------------------------------


# Load the correct marker database
load(file.path(geno_dir, "phased_marker_genotype_db.RData"))

# Data.frame of snp metadata
snp_info <- phased_geno_hmp %>%
  select(marker, chrom, pos, alleles)

# Vectors of germplasm from the three different categories in pop_metadata
germplasm_names_list <- pop_metadata %>%
  split(.$category) %>%
  map("individual") %>%
  map(~intersect(., names(phased_geno_hmp)))

# Subset the marker data for the wild germplasm
geno_mat_wild <- phased_geno_mat[germplasm_names_list$Wild,] - 1
haplo_array_wild <- phased_geno_haplo_array[,,germplasm_names_list$Wild]

geno_mat_native <- phased_geno_mat[germplasm_names_list$NativeSelection,] - 1
haplo_array_native <- phased_geno_haplo_array[,,germplasm_names_list$NativeSelection]

geno_mat_breeding <- phased_geno_mat[germplasm_names_list$Breeding,] - 1
haplo_array_breeding <- phased_geno_haplo_array[,,germplasm_names_list$Breeding]

# Calculate MAF for each group
maf_wild <- calc_maf(x = geno_mat_wild)
maf_native <- calc_maf(x = geno_mat_native)
maf_breeding <- calc_maf(x = geno_mat_breeding)

# Plot together
par(mfrow = c(2,2))
hist(maf_wild, main = "MAF wild accessions - prefilter")
hist(maf_native, main = "MAF native accessions - prefilter")
hist(maf_breeding, main = "MAF breeding accessions - prefilter")
par(mfrow = c(1,1))


# Remove identical (i.e. all hets)
genotypes_per_marker <- apply(X = geno_mat_wild, MARGIN = 2, FUN = n_distinct)
geno_mat_wild1 <- geno_mat_wild[,genotypes_per_marker > 1, drop = FALSE]
maf_wild <- calc_maf(x = geno_mat_wild1)


# First filter out SNPs that are monomorphic or in complete LD
geno_mat_wild_filter1 <- prune_LD(x = geno_mat_wild1[, maf_wild > 0], r2.max = max_LD_r2)
geno_mat_wild_filter1 <- geno_mat_wild_filter1[,!is.na(colnames(geno_mat_wild_filter1))]

# Dimensions
dim(geno_mat_wild_filter1)

## Calculate relationship matrices
# Wild cranberry
K_wild <- A.mat(X = geno_mat_wild_filter1, min.MAF = 0, max.missing = 1)

# Calculate the K matrix with wild and native selection samples
K_all <- A.mat(X = rbind(geno_mat_wild_filter1, rbind(geno_mat_native, geno_mat_breeding)[,colnames(geno_mat_wild_filter1)]),
               min.MAF = 0, max.missing = 1)


# Next filter on MAF - this will be used to calculate LD
geno_mat_wild_filter2 <- filter_snps(x = geno_mat_wild_filter1, r2.max = max_LD_r2, maf.min = min_mac / nrow(geno_mat_wild_filter1),
                                     indiv.miss.max = max_indiv_miss, snp.miss.max = max_snp_miss)
geno_mat_wild_filter2 <- geno_mat_wild_filter2[,!is.na(colnames(geno_mat_wild_filter2))]

# Dimensions
dim(geno_mat_wild_filter2)

# Now filter on MAF AND a lower LD threshold - these will be used for mapping
geno_mat_wild_filter3 <- filter_snps(x = geno_mat_wild_filter1, r2.max = max_LD_r2_2,
                                     maf.min = min_mac / nrow(geno_mat_wild_filter1),
                                     indiv.miss.max = max_indiv_miss, snp.miss.max = max_snp_miss)
geno_mat_wild_filter3 <- geno_mat_wild_filter3[,!is.na(colnames(geno_mat_wild_filter3))]

# Dimensions
dim(geno_mat_wild_filter3)




# Recalculate and visualize MAF
maf <- calc_maf(x = geno_mat_wild_filter3); hist(maf, main = "MAF - postfilter")

# Filter for those same SNPs in the haplo array
haplo_array_wild_filter3 <- haplo_array_wild[,colnames(geno_mat_wild_filter3),]

## Filter the entire marker genotype matrix and haplotype matrix for these markers
geno_mat_all_filter3 <- phased_geno_mat[row.names(K_all),colnames(geno_mat_wild_filter3)] - 1
haplo_array_all_filter3 <- phased_geno_haplo_array[,colnames(geno_mat_wild_filter3),row.names(K_all)]


# Rename and save
geno_hmp_wild_filter3 <- t(geno_mat_wild_filter3) %>%
  as.data.frame() %>%
  rownames_to_column("marker") %>%
  filter(marker %in% colnames(geno_mat_wild_filter3)) %>%
  inner_join(snp_info, .) %>%
  select(marker:alleles, row.names(geno_mat_wild_filter3))


# Filter snp_info
snp_info_all <- snp_info
snp_info <- snp_info_all %>%
  filter(marker %in% colnames(geno_mat_wild_filter3))


save("geno_mat_wild_filter2", "geno_mat_wild_filter3", "haplo_array_wild_filter3",
     "geno_mat_all_filter3", "haplo_array_all_filter3", "geno_hmp_wild_filter3",
     "K_wild", "K_all", "snp_info", "snp_info_all", "pop_metadata",
     file = file.path(data_dir, "population_metadata_and_genotypes.RData"))




## Export the genotype data as a vcf
data("vcfR_test")

# First create a data.frame of the fixed section
vcf_fixed <- snp_info %>%
  filter(marker %in% colnames(haplo_array_wild_filter3)) %>%
  separate(alleles, c("ref", "alt"), sep = "/") %>%
  rename_all(toupper) %>%
  select(CHROM, POS, ID = MARKER, names(.)) %>%
  mutate(QUAL = ".", FILTER = "PASS", INFO = ".") %>%
  mutate_all(as.character) %>%
  as.matrix()

# Create the genotype section
vcf_gt <- apply(X = haplo_array_wild_filter3, MARGIN = c(2, 3), FUN = function(genotype) paste0(genotype, collapse = "|") ) %>%
  as.data.frame() %>%
  .[vcf_fixed[,"ID"],] %>%
  cbind(FORMAT = "GT", .) %>%
  as.matrix()

# Create the meta vector
vcf_meta <- c(
  "##fileformat=VCFv4.0",
  paste0("##fileDate=", format(Sys.Date(), "%Y%m%d")),
  "##source=R",
  "##phasing=full",
  "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
)

# Assemble the VCF
vcf_wild <- vcfR_test
vcf_wild@meta <- vcf_meta
vcf_wild@fix <- vcf_fixed
vcf_wild@gt <- vcf_gt

# Write the VCF
vcf_file <- file.path(data_dir, "wild_cranberry_cleaned_genotypes_filter3.vcf.gz")
write.vcf(x = vcf_wild, file = vcf_file)


## Convert the VCF to PLINK BED format
# First read the vcf using snprelate
out_gds_file <- "temp.gds"
gds <- SNPRelate::snpgdsVCF2GDS(vcf.fn = vcf_file, out.fn = out_gds_file)
# Read in the gds file
gds_in <- SNPRelate::snpgdsOpen(filename = gds)
# Convert to plink
SNPRelate::snpgdsGDS2PED(gdsobj = gds_in, ped.fn = str_remove(vcf_file, ".vcf.gz"))
# Delete the temp file
SNPRelate::snpgdsClose(gds_in)
file.remove(gds)

