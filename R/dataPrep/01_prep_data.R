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

# Directory of phenotypic data
pheno_dir <- file.path(cran_dir, "Breeding/prior2021/PhenotypeData/PhenotypicAnalysis/results/")
# Directory of genotypic data
geno_dir <- file.path(cran_dir, "Genotyping/MarkerDatabase")


## User parameters
max_LD_r2 <- 0.99999 # Max LD of any one pair of markers
min_mac <- 10 # Minimum minor allele count
max_indiv_miss <- 0.7
max_snp_miss <- 0.7


# Copy phenotypic data ----------------------------------------------------


# Read in the germplasm collection phenotypic data
filename <- file.path(pheno_dir, "germplasm_collection_phenotypic_blups.xlsx")
pheno_dat <- read_excel(path = filename, sheet = "blups")


# Copy germplasm collection metadata --------------------------------------

# Path to the filename
filename <- file.path(cran_dir, "Breeding/prior2021/Populations/GermplasmCollection/germplasmCollectionMetadata.xlsx")
metadata <- read_excel(path = filename, na = c("", "NA"))

# Read in variety descriptions
filename <- file.path(cran_dir, "CranberryInformation/VarietyDescriptions/cranberry_variety_descriptions_quantified.xlsx")
variety_desc <- read_excel(path = filename, na = c("", "NA")) %>%
  mutate(possible_germplasm_collection_id = str_split(possible_germplasm_collection_id, ", "))


# Subset those with marker data; select relevant columns
metadata1 <- metadata %>%
  filter(!is.na(marker_sample_name)) %>%
  mutate(genetic_id = str_split(genetic_id, ", ")) %>%
  select(marker_sample_name, individual, formatted_name, genetic_id, variety_designation, variety_selection_name,
         category, contains("origin"))

# Vector of wild germplasm accession names
wild_germplasm_names <- subset(metadata1, category == "wild", individual, drop = TRUE)
# Vector of native selection germplasm accession names
native_selection_germplasm_names <- subset(metadata1, category == "native_selection", individual, drop = TRUE)
# Vector of cultivar germplasm accession names
cultivar_germplasm_names <- subset(metadata1, category == "cultivar_gen1", individual, drop = TRUE)


## Save these to an RData file
germ_meta <- metadata1 %>%
  mutate(genetic_id = map_chr(genetic_id, ~paste0(., collapse = ", ")) %>% parse_character())


save("pheno_dat", "germ_meta", "wild_germplasm_names", "native_selection_germplasm_names",
     "cultivar_germplasm_names", file = file.path(data_dir, "germplasm_metadata.RData"))


# Copy marker genotype data -----------------------------------------------


# Load the correct marker database
load(file.path(geno_dir, "phased_marker_genotype_db.RData"))

# Data.frame of snp metadata
snp_info <- phased_geno_hmp %>%
  select(marker, chrom, pos, alleles)

# Intersect wild germplasm names with names in the marker matrix
wild_germplasm_names_genotyped <- intersect(wild_germplasm_names, row.names(phased_geno_mat))
# Same with native selections
native_selection_germplasm_names_genotyped <- intersect(native_selection_germplasm_names, row.names(phased_geno_mat))

## Filter the variety description df for germplasm accessions that were genotyped
## Save this new df
variety_desc1 <- variety_desc %>%
  filter(map_lgl(possible_germplasm_collection_id, ~any(. %in% native_selection_germplasm_names_genotyped))) %>%
  mutate(possible_germplasm_collection_id = map(possible_germplasm_collection_id, intersect, native_selection_germplasm_names_genotyped)) %>%
  select(variety_name = variety, germplasm_collection_name = possible_germplasm_collection_id, original_planting_state) %>%
  unnest(germplasm_collection_name) %>%
  write_csv(x = ., file = file.path(data_dir, "native_selection_metadata.csv"))





# Subset the marker data for the wild germplasm
geno_mat_wild <- phased_geno_mat[wild_germplasm_names_genotyped,] - 1
haplo_array_wild <- phased_geno_haplo_array[,,wild_germplasm_names_genotyped]

# First filter for LD and monomorphic SNPs
maf <- calc_maf(x = geno_mat_wild); hist(maf, main = "MAF - prefilter")
geno_mat_wild1 <- prune_LD(x = geno_mat_wild[, maf > 0], r2.max = max_LD_r2)
geno_mat_wild1 <- geno_mat_wild1[,!is.na(colnames(geno_mat_wild1))]


## Calculate relationship matrices
# Wild cranberry
K_wild <- A.mat(X = geno_mat_wild1, min.MAF = 0, max.missing = 1)
# All genotypes - leave-one-chromosome-out
K_wild_loco <- snp_info %>%
  filter(marker %in% colnames(geno_mat_wild1)) %>%
  split(.$chrom) %>%
  map("marker") %>%
  map(~setdiff(colnames(geno_mat_wild1), .x)) %>%
  map(~A.mat(X = geno_mat_wild1[, .x, drop = FALSE],
             min.MAF = 0, max.missing = 1))

# Calculate the K matrix with wild and native selection samples
K_all <- A.mat(X = rbind(geno_mat_wild1, phased_geno_mat[native_selection_germplasm_names_genotyped, colnames(geno_mat_wild1)] - 1),
               min.MAF = 0, max.missing = 1)


# Next filter those SNPs for higher MAF
geno_mat_wild2 <- filter_snps(x = geno_mat_wild1, r2.max = max_LD_r2, maf.min = min_mac / nrow(geno_mat_wild1),
                              indiv.miss.max = max_indiv_miss, snp.miss.max = max_snp_miss)

# Remove identical (i.e. all hets)
genotypes_per_marker <- apply(X = geno_mat_wild2, MARGIN = 2, FUN = n_distinct)
geno_mat_wild2 <- geno_mat_wild2[,genotypes_per_marker > 1, drop = FALSE]


# Recalculate and visualize MAF
maf <- calc_maf(x = geno_mat_wild2); hist(maf, main = "MAF - postfilter")

# Filter for those same SNPs in the haplo array
haplo_array_wild1 <- haplo_array_wild[,colnames(geno_mat_wild2),]

## Filter the entire marker genotype matrix and haplotype matrix for these markers
geno_mat_all <- phased_geno_mat[,colnames(geno_mat_wild2)] - 1
haplo_array_all <- phased_geno_haplo_array[,colnames(geno_mat_wild2),]


# Rename and save
geno_mat_wild <- geno_mat_wild2
haplo_array_wild <- haplo_array_wild1
geno_hmp_wild <- t(geno_mat_wild) %>%
  as.data.frame() %>%
  rownames_to_column("marker") %>%
  filter(marker %in% colnames(geno_mat_wild)) %>%
  inner_join(snp_info, .) %>%
  select(marker:alleles, row.names(geno_mat_wild))

# Keep a version with no filtering
geno_mat_unfiltered <- phased_geno_mat[c(wild_germplasm_names_genotyped, native_selection_germplasm_names_genotyped),] - 1



save("geno_mat_wild", "haplo_array_wild", "geno_mat_all", "haplo_array_all", "geno_hmp_wild", "geno_mat_unfiltered",
     "K_wild", "K_wild_loco", "K_all", "snp_info",
     file = file.path(data_dir, "gc_marker_data.RData"))




## Export the genotype data as a vcf
data("vcfR_test")

# First create a data.frame of the fixed section
vcf_fixed <- snp_info %>%
  filter(marker %in% colnames(haplo_array_wild)) %>%
  separate(alleles, c("ref", "alt"), sep = "/") %>%
  rename_all(toupper) %>%
  select(CHROM, POS, ID = MARKER, names(.)) %>%
  mutate(QUAL = ".", FILTER = "PASS", INFO = ".") %>%
  mutate_all(as.character) %>%
  as.matrix()

# Create the genotype section
vcf_gt <- apply(X = haplo_array_wild, MARGIN = c(2, 3), FUN = function(genotype) paste0(genotype, collapse = "|") ) %>%
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
vcf_file <- file.path(data_dir, "wild_cranberry_cleaned_genotypes.vcf.gz")
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

