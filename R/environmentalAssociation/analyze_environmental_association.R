# Germplasm collection environmental association
#
# Analyze environmental association results
#

# Load packages
library(sommer)
library(raster)
library(neyhart)
library(tidyverse)
library(readxl)
library(snps)
library(qvalue)
library(LDheatmap)
library(slider)
library(patchwork)
library(cowplot)

# Load the startup script
source("startup.R")


# Read in the marker data
load(file.path(data_dir, "gc_marker_data.RData"))
load(file.path(data_dir, "germplasm_metadata.RData"))

# Load the worldclim data
load(file.path(data_dir, "germplasm_origin_bioclim_data.RData"))

# Center and scale
eaa_environmental_data1 <- eaa_environmental_data %>%
  select(location_abbr, all_of(eaa_environmental_vars$variable)) %>%
  # mutate_at(vars(contains("bio")), scale, scale = FALSE) %>%
  # mutate_at(vars(contains("bio")), as.numeric) %>%
  mutate(elevation = ifelse(is.na(elevation), 0, elevation))

# Read in the germplasm metadata
germ_meta <- germ_meta %>%
  left_join(., select(eaa_environmental_data, origin_name, location_abbr)) %>%
  # Get the state/province name from the location
  mutate(state = str_sub(string = location_abbr, start = -2, end = -1))


# Read in variety origin information
native_sel_meta <- read_csv(file = file.path(data_dir, "native_selection_metadata.csv")) %>%
  rename(state = original_planting_state)


# Load the eGWAS results
load(file.path(result_dir, "eaa_gwas_results.RData"))
# Load the pop gen results
load(file.path(result_dir, "population_genetics_stats.RData"))



# Prepare the bioclim data ------------------------------------------------

# Prepare the bioclim data as if it were trait data
germplasm_bioclim_data <- germ_meta %>%
  filter(category == "wild", !is.na(origin_latitude)) %>%
  select(individual, latitude = origin_latitude, longitude = origin_longitude, state) %>%
  left_join(., eaa_environmental_data1) %>%
  as.data.frame() %>%
  filter(individual %in% colnames(K_wild))




# Manipulate the marker data -------------------------------------------------

geno_mat_wild1 <- geno_mat_wild[germplasm_bioclim_data$individual,]

geno_hmp_wild1 <- geno_hmp_wild %>%
  select(marker:alleles, row.names(geno_mat_wild1))

snp_info1 <- snp_info %>%
  rename(Locus = marker, LG = chrom, Position = pos) %>%
  as.data.frame()

# Split the marker matrix by chromosome
geno_mat_list <- geno_hmp_wild %>%
  select(marker, chrom) %>%
  split(.$chrom) %>% map("marker") %>%
  map(~geno_mat_wild1[,.x, drop = FALSE])

# Subset the geno_hmp_wild object
geno <- geno_hmp_wild %>%
  select(marker, chrom, pos, all_of(colnames(K_wild))) %>%
  mutate_at(vars(colnames(K_wild)), ~. - 1) %>%
  as.data.frame()


# Subset the native selections
# geno_mat_native <- geno_mat_all[native_sel_meta$germplasm_collection_name,]
geno_mat_native <- geno_mat_all[intersect(row.names(geno_mat_all), native_selection_germplasm_names),]
geno_mat_cultivar <- geno_mat_all[intersect(row.names(geno_mat_all), cultivar_germplasm_names),]






# Compare significant marker alleles with bioclim variables ---------------

# Extract the significant markers from the multi-locus mixed model output
egwas_sigmar <- eaa_gwas_best_model_signif_marker %>%
  unnest(sigmar) %>%
  mutate(class = "egwas_sigmar") %>% filter(! variable %in% c("elevation", "latitude", "longitude", "IC1", "IC2", "IC3"))

# Number of unique markers per variable
egwas_sigmar %>%
  left_join(., select(eaa_environmental_vars, -class)) %>%
  xtabs(~ full_name, .) %>%
  sort()

# Number of unique markers
egwas_sigmar_unique <- egwas_sigmar %>%
  select(marker, variable) %>%
  group_by(marker) %>%
  nest(data = variable) %>%
  mutate(variable = map(data, "variable")) %>%
  arrange(desc(map_dbl(variable, length))) %>%
  mutate(nVar = map_dbl(data, nrow),
         variable = map_chr(variable, ~paste0(., collapse = ", "))) %>%
  ungroup() %>%
  select(marker, nVar, variable)


# Tranpose the genotype hmp object
geno_hmp_wild2 <- geno_hmp_wild1 %>%
  separate(col = alleles, c("ref_allele", "alt_allele"), sep = "/") %>%
  mutate(marker1 = paste0(marker, "_", ref_allele, "_", alt_allele)) %>%
  filter(marker %in% unique(egwas_sigmar$marker)) %>%
  mutate_at(vars(-contains("marker"), -marker:-alt_allele), ~. + 1) %>%
  mutate_at(vars(-contains("marker"), -marker:-alt_allele), ~pmap_chr(list(., ref_allele, alt_allele), ~{
    snp <- ..1
    ref <- ..2
    alt <- ..3
    paste0(rep(ref, snp), rep(alt, 2 - snp), collapse = "")

  })) %>%
  select(contains("marker"), chrom:alt_allele, names(.)) %>%
  gather(individual, genotype, -marker:-alt_allele)


geno_hmp_wild1t <- geno_hmp_wild2 %>%
  select(marker, individual, genotype) %>%
  spread(marker, genotype)

geno_hmp_wild1t2 <- geno_hmp_wild1 %>%
  mutate_at(vars(-contains("marker"), -marker:-alleles), ~. + 1) %>%
  gather(individual, genotype, -marker:-alleles) %>%
  select(marker, individual, genotype) %>%
  spread(marker, genotype)

# Marker metadata
marker_metadata <- geno_hmp_wild2 %>%
  distinct_at(vars(marker, chrom, pos, contains("allele")))


# Allele colors
allele_colors <- setNames(neyhart_palette("umn2")[3:4], c("ref_allele", "alt_allele"))



# Compare eGWAS results with SPA --------------------------------------

# Read in the SPA results
spa_results <- spa_out


# Distribution of spa scores
spa_results %>%
  ggplot(aes(x = spa_score)) +
  geom_histogram(bins = 50) +
  theme_classic()

# Get the 1% quantile of spa scores
spa_scores_q001 <- quantile(spa_results$spa_score, 1 - 0.01)


# Get the spa_scores for the significant GWAS hits
eaa_gwas_spa_scores <- eaa_gwas_best_model_signif_marker %>%
  unnest(sigmar) %>%
  distinct(marker) %>%
  left_join(., select(spa_results, marker, spa_score))

# How many of these are above the threshold?
mean(eaa_gwas_spa_scores$spa_score >= spa_scores_q001)

# Arrange in descending spa score order
eaa_gwas_spa_scores %>%
  arrange(desc(spa_score)) %>%
  left_join(., egwas_sigmar_unique)

eaa_gwas_spa_scores %>%
  arrange(desc(spa_score)) %>%
  left_join(., egwas_sigmar_unique) %>%
  arrange(desc(nVar))

# Average spa score of GWAS hits
gwas_mean_spa <- mean(eaa_gwas_spa_scores$spa_score)
hist(eaa_gwas_spa_scores$spa_score)


## Randomly sample the same number of background SNPs as GWAS SNPs and compare
## mean SPA scores
spa_out_nogwas_vec <- spa_out %>%
  subset(x = ., ! marker %in% eaa_gwas_spa_scores$marker, spa_score, drop = TRUE)

random_snps_spa <- replicate(n = 10000, mean(sample(x = spa_out_nogwas_vec, size = nrow(eaa_gwas_spa_scores))))

# P-value for the GWAS SNPs SPA score
mean(random_snps_spa >= gwas_mean_spa)




# Look for gene annotations ------------------------------------------------

# read in the GFF file
cranberry_gff <- ape::read.gff(file = file.path(data_dir, "Vaccinium_macrocarpon_BenLear_v2_annotations.gff"))


## Look for gene annotations near the significant markers

# What bp distance should we use?
bp <- 17000 # Use 17000; this is the bp at which LD decays to 0.2

# Iterate over the signficant markers
egwas_sigmar_nearby_annotation <- egwas_sigmar %>%
  group_by(variable, marker) %>%
  do(nearby_annotation = {
    row <- .

    # Get the marker location
    snp_info_i <- subset(snp_info, marker == row$marker)

    # Filter the GFF for annotations near the marker
    sigmar_nearby_ann <- cranberry_gff %>%
      filter(seqid == as.numeric(snp_info_i$chrom)) %>%
      filter(start >= snp_info_i$pos - bp,
             end <= snp_info_i$pos + bp) %>%
      as_tibble() %>%
      mutate_at(vars(start, end), list(distance = ~abs(. - snp_info_i$pos))) %>%
      # Parse out the attributes
      mutate(attributes = str_split(attributes, ";"),
             attributes = map(attributes, ~{
               splt <- str_split(.x, "=")
               as_tibble(setNames(map(splt, 2), map_chr(splt, 1)))
             })) %>%
      select(type, contains("start"), contains("end"), score:attributes) %>%
      arrange(start_distance)

    sigmar_nearby_ann

  }) %>% ungroup()

# How many nearby genes?
egwas_sigmar_nearby_annotation1 <- egwas_sigmar_nearby_annotation %>%
  mutate(nGenes = map_dbl(nearby_annotation, ~nrow(subset(.x, type == "gene"))),
         closest_gene_dist = map_dbl(nearby_annotation, ~subset(.x, type == "gene", start_distance, drop = T)[1]))



# How close is each SNP marker to its nearest gene?
all_marker_gene_distance <- snp_info %>%
  filter(marker %in% colnames(geno_mat_wild)) %>%
  mutate(chrom = parse_number(chrom)) %>%
  # head(250) %>%
  group_by(marker) %>%
  do({
    row <- .
    subset(cranberry_gff, seqid == row$chrom & type == "gene") %>%
      mutate_at(vars(start, end), list(distance = ~abs(.x - row$pos))) %>%
      mutate(min_distance = pmin(start_distance, end_distance)) %>%
      top_n(x = ., n = 1, wt = -min_distance) %>%
      select(min_distance)

  }) %>% ungroup()


egwas_sigmar_nearby_annotation1


# Save the annotation data
save("egwas_sigmar_nearby_annotation1", "egwas_sigmar",
     file = file.path(result_dir, "egwas_sigmar_analysis.RData"))























