# CranberryGermplasmEAA
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

# FDR threshold
fdr_thresh <- 0.20


# Read in the base data
load(file.path(data_dir, "population_metadata_and_genotypes.RData"))
# Load the worldclim data
load(file.path(data_dir, "germplasm_origin_bioclim_data.RData"))

# Load the eGWAS results
load(file.path(result_dir, "eaa_gwas_results.RData"))
# Load the pop gen results
load(file.path(result_dir, "population_genetics_stats.RData"))


# Center and scale
eaa_environmental_data1 <- eaa_environmental_data %>%
  select(location_abbr, all_of(eaa_environmental_vars$variable))


# Rename genotype matrices
geno_mat_wild <- geno_mat_wild_filter3
geno_hmp_wild <- geno_hmp_wild_filter3
geno_mat_all <- geno_mat_all_filter3


# Prepare the bioclim data ------------------------------------------------

# Prepare the bioclim data as if it were trait data
germplasm_bioclim_data <- pop_metadata %>%
  filter(category == "Wild", !is.na(latitude)) %>%
  select(individual, latitude, longitude) %>%
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




# Compare significant marker alleles with bioclim variables ---------------

# Subset GWAS results for the best model
gwas_out_best_model <- gwas_out %>%
  inner_join(., subset(p_lambda, best_model == "*", c(model, variable), drop = TRUE))

# Count the best models
xtabs(~ model, p_lambda, subset = best_model == "*")

# Use model2
gwas_out_best_model <- gwas_out %>% filter(model == "model2")

# Extract the significant markers based on FDR
egwas_sigmar <- eaa_gwas_sigmar

# Number of total associations
nrow(egwas_sigmar)

# 126

# Number of unique markers
n_distinct(egwas_sigmar$marker)

# 60

n_distinct(egwas_sigmar$variable)

# 34


# Number of unique markers per variable
egwas_sigmar %>%
  left_join(., select(eaa_environmental_vars, -class)) %>%
  subset(variable %in% env_variables_of_interest) %>%
  xtabs(~ full_name, .) %>%
  sort()

# 1 - 9

# Number of variables per marker
egwas_sigmar %>%
  left_join(., select(eaa_environmental_vars, -class)) %>%
  subset(variable %in% env_variables_of_interest) %>%
  xtabs(~ marker, .) %>%
  sort()

# 1 - 8

# Number of markers per variable class
egwas_sigmar_cat <- egwas_sigmar %>%
  left_join(., eaa_environmental_vars) %>%
  subset(variable %in% env_variables_of_interest)

egwas_sigmar_cat %>%
  group_by(class, variable) %>%
  summarize(nSigMar = n()) %>%
  summarize(nVariable = n_distinct(variable), nSigMar = sum(nSigMar))

# class               nVariable nSigMar
# 1 geography                   3       7
# 2 precipitation               7      25
# 3 principal_component         2       7
# 4 soil                       15      68
# 5 temperature                 7      19


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


# Get the 1% quantile of spa scores
spa_scores_q001 <- quantile(spa_results$spa_score, 0.999)

# Subset markers with spa scores above this threshold
spa_outliers <- spa_results %>%
  filter(spa_score >= spa_scores_q001)

# Same with Fst results
fst_q001 <- quantile(global_marker_fst$Fst, 0.999)
fst_outliers <- global_marker_fst %>%
  filter(Fst >= fst_q001)

# Overlap between Fst and SPA
intersect(fst_outliers$marker, spa_outliers$marker)

# None


# Get the spa_scores for the significant GWAS hits
eaa_gwas_spa_scores <- egwas_sigmar %>%
  distinct(marker) %>%
  left_join(., select(spa_results, marker, spa_score))

# How many of these are above the threshold?
sum(eaa_gwas_spa_scores$spa_score >= spa_scores_q001)

# Just 1
#
subset(eaa_gwas_spa_scores, spa_score >= spa_scores_q001)

# Arrange in descending spa score order
eaa_gwas_spa_scores %>%
  arrange(desc(spa_score)) %>%
  left_join(., egwas_sigmar) %>%
  # Calculate the percentile of each marker spa
  mutate(spa_percentile = map_dbl(spa_score, ~mean(. >= spa_results$spa_score)))

# Average spa score of GWAS hits
gwas_mean_spa <- mean(eaa_gwas_spa_scores$spa_score)
sd(eaa_gwas_spa_scores$spa_score)

# Mean and SD of all marker SPA scores
mean(spa_out$spa_score)
sd(spa_out$spa_score)

hist(eaa_gwas_spa_scores$spa_score)
# Number of unique egwas markers
n_eaa_markers <- length(eaa_gwas_spa_scores$spa_score)


## Randomly sample the same number of background SNPs as GWAS SNPs and compare
## mean SPA scores
spa_out_nogwas_vec <- subset(x = spa_out, ! marker %in% eaa_gwas_spa_scores$marker, spa_score, drop = TRUE)



set.seed(209)
random_snps_spa <- replicate(n = 10000, mean(sample(x = spa_out_nogwas_vec, size = n_eaa_markers)))

# P-value for the GWAS SNPs SPA score
mean(random_snps_spa >= gwas_mean_spa)
hist(random_snps_spa, breaks = 50); abline(v = gwas_mean_spa, col = "blue")

wilcox.test(x = eaa_gwas_spa_scores$spa_score, y = spa_out_nogwas_vec)

# P = 0.01345

# Save this as a figure
g_random_eaa_spa <- ggplot(data = NULL, aes(x = random_snps_spa)) +
  geom_histogram(bins = 50, fill = "gray70") +
  geom_vline(xintercept = gwas_mean_spa, color = "blue") +
  geom_text(aes(x = gwas_mean_spa, y = 600, label = "EAA signif.\nSNPs"), color = "blue",
            nudge_x = 0.02, hjust = 0, size = 3) +
  geom_text(aes(x = min(random_snps_spa), y = 600, label = paste0("Random sets of\n", n_eaa_markers, " SNPS")),
            color = "gray70", hjust = 1, nudge_x = 0.1, size = 3) +
  scale_x_continuous(name = "Mean SPA score", breaks = pretty) +
  scale_y_continuous(name = "Frequency", breaks = pretty) +
  theme_genetics()

ggsave(filename = "FigureSXX_eaaSPA_versus_randomSPA.jpg", plot = g_random_eaa_spa, path = fig_dir,
       height = 2, width = 4, dpi = 500)

# Correlation between eaa p-values and SPA
egwas_score_spa <- gwas_out_best_model %>%
  filter(variable %in% unique(egwas_sigmar$variable)) %>%
  left_join(., select(spa_out, marker, spa_score)) %>%
  group_by(variable) %>%
  summarize(spa_egwas_cor = cor(sqrt(score), spa_score), .groups = "drop")




### Do the same thing with fst ###

# Get the fst for the significant GWAS hits
eaa_gwas_fst <- eaa_gwas_sigmar %>%
  distinct(marker) %>%
  left_join(., select(global_marker_fst, marker, Fst))

# How many of these are above the threshold?
sum(eaa_gwas_fst$Fst >= fst_q001)

# 2

# Number of unique egwas markers
n_eaa_markers <- length(eaa_gwas_fst$Fst)


# 60

# Intersect eaa_gwas markers and Fst outliers
intersect(fst_outliers$marker, unique(eaa_gwas_sigmar$marker))

# 2 markers overlap

# Arrange in descending Fst order
eaa_gwas_fst %>%
  arrange(desc(Fst)) %>%
  left_join(., egwas_sigmar) %>%
  # Calculate the percentile of each marker Fst
  mutate(fst_percentile = map_dbl(Fst, ~mean(. >= global_marker_fst$Fst)))

# Average spa score of GWAS hits
gwas_mean_fst <- mean(eaa_gwas_fst$Fst)
hist(eaa_gwas_fst$Fst)

## Randomly sample the same number of background SNPs as GWAS SNPs and compare
## mean SPA scores
fst_nogwas_vec <- subset(x = global_marker_fst, ! marker %in% eaa_gwas_fst$marker, Fst, drop = TRUE)

set.seed(209)
random_snps_fst <- replicate(n = 10000, mean(sample(x = fst_nogwas_vec, size = nrow(eaa_gwas_fst))))

# P-value for the GWAS SNPs SPA score
mean(random_snps_fst >= gwas_mean_fst)
# essentially p = 0

hist(random_snps_fst, breaks = 50, xlim = range(pretty(c(random_snps_fst, gwas_mean_fst)))); abline(v = gwas_mean_fst, col = "blue")

wilcox.test(x = eaa_gwas_fst$Fst, y = fst_nogwas_vec)

# p-value < 2.2e-16

# Save this as a figure
g_random_eaa_fst <- ggplot(data = NULL, aes(x = random_snps_fst)) +
  geom_histogram(bins = 50, fill = "gray70") +
  # geom_histogram(aes(x = eaa_gwas_fst$Fst), fill = "red") +
  geom_vline(xintercept = gwas_mean_fst, color = "red") +
  geom_text(aes(x = gwas_mean_fst, y = 1000, label = "EAA signif.\nSNPs"), color = "red",
            nudge_x = -0.02, hjust = 1, size = 3) +
  geom_text(aes(x = min(random_snps_fst), y = 1000, label = paste0("Random sets of\n", n_eaa_markers, " SNPS")),
            color = "gray70", hjust = 0, nudge_x = 0.08, size = 3) +
  scale_x_continuous(name = expression('Mean'~italic(F)[ST]), breaks = pretty) +
  scale_y_continuous(name = "Frequency", breaks = pretty) +
  theme_genetics()

ggsave(filename = "FigureSXX_eaaFST_versus_randomFST.jpg", plot = g_random_eaa_fst, path = fig_dir,
       height = 2, width = 4, dpi = 500)


# Correlation between eaa p-values and Fst
egwas_score_fst <- gwas_out_best_model %>%
  filter(variable %in% unique(egwas_sigmar$variable)) %>%
  left_join(., select(global_marker_fst, marker, Fst)) %>%
  group_by(variable) %>%
  summarize(fst_egwas_cor = cor(sqrt(score), Fst), .groups = "drop")

egwas_score_fst %>%
  arrange(desc(abs(fst_egwas_cor)))







# Look for gene annotations ------------------------------------------------

gff_file <- file.path(data_dir, "/Vaccinium_macrocarpon_BenLear_v2_annotations.gff")

# read in the GFF file
cranberry_gff <- ape::read.gff(file = gff_file)

## Look for gene annotations near the significant markers

# What bp distance should we use?
bp <- 17000 # Use 17000; this is the bp at which LD decays to 0.2

# Combine egwas, Fst, and SPA markers
all_sigmars <- bind_rows(
  mutate(egwas_sigmar, class = "EAA"),
  mutate(fst_outliers, variable = "Fst", class = variable),
  mutate(spa_outliers, variable = "SPA", class = variable)
) %>% select(class, variable, marker, chrom, pos)

# Iterate over the significant markers
all_sigmar_nearby_annotation <- all_sigmars %>%
  group_by(class, variable, marker) %>%
  do({
    row <- .

    # Get the marker location
    snp_info_i <- subset(snp_info, marker == row$marker)

    # Is this SNP within any genes?
    snp_i_within_annotation <- cranberry_gff %>%
      filter(seqid == as.numeric(snp_info_i$chrom),
             type != "chromosome") %>%
      filter((start <= snp_info_i$pos & end >= snp_info_i$pos) | (start >= snp_info_i$pos & end <= snp_info_i$pos))

    # Filter the GFF for annotations near the marker
    sigmar_nearby_ann <- cranberry_gff %>%
      filter(seqid == as.numeric(snp_info_i$chrom), type != "chromosome") %>%
      filter(abs(snp_info_i$pos - start) <= bp | abs(snp_info_i$pos - end) <= bp) %>%
      as_tibble() %>%
      mutate_at(vars(start, end), list(distance = ~abs(. - snp_info_i$pos))) %>%
      # Parse out the attributes
      mutate(attributes = str_split(attributes, ";"),
             attributes = map(attributes, ~{
               splt <- str_split(.x, "=")
               as_tibble(setNames(map(splt, 2), map_chr(splt, 1)))
             }),
             min_distance = pmin(start_distance, end_distance)) %>%
      select(type, contains("start"), contains("end"), score:attributes, min_distance) %>%
      arrange(min_distance)

    tibble(nearby_annotation = list(sigmar_nearby_ann), within_annotations = list(snp_i_within_annotation))

  }) %>% ungroup() %>%
  mutate(region = ifelse(sapply(within_annotations, nrow) > 0, "genic", "non-genic"))


# Iterate over all markers
all_markers_nearby_annotation <- snp_info %>%
  group_by(marker) %>%
  do({
    snp_info_i <- .

    # Is this SNP within any genes?
    snp_i_within_annotation <- cranberry_gff %>%
      filter(seqid == as.numeric(snp_info_i$chrom),
             type != "chromosome") %>%
      filter((start <= snp_info_i$pos & end >= snp_info_i$pos) | (start >= snp_info_i$pos & end <= snp_info_i$pos))

    # Filter the GFF for annotations near the marker
    sigmar_nearby_ann <- cranberry_gff %>%
      filter(seqid == as.numeric(snp_info_i$chrom), type != "chromosome") %>%
      filter(abs(snp_info_i$pos - start) <= bp | abs(snp_info_i$pos - end) <= bp) %>%
      as_tibble() %>%
      mutate_at(vars(start, end), list(distance = ~abs(. - snp_info_i$pos))) %>%
      # Parse out the attributes
      mutate(attributes = str_split(attributes, ";"),
             attributes = map(attributes, ~{
               splt <- str_split(.x, "=")
               as_tibble(setNames(map(splt, 2), map_chr(splt, 1)))
             }),
             min_distance = pmin(start_distance, end_distance)) %>%
      select(type, contains("start"), contains("end"), score:attributes, min_distance) %>%
      arrange(min_distance)

    tibble(min_distance_nearby_gene = min(subset(sigmar_nearby_ann, type == "gene", min_distance, drop = TRUE)),
           region = ifelse(nrow(snp_i_within_annotation) == 0, "non-genic", "genic"))

  }) %>% ungroup()


# Replace Inf with NA
all_markers_nearby_annotation <- all_markers_nearby_annotation %>%
  mutate(min_distance_nearby_gene = ifelse(is.infinite(min_distance_nearby_gene), NA, min_distance_nearby_gene))


# Proportion of genic SNPs in significant SNPs versus all SNPs
all_sigmar_nearby_annotation %>%
  distinct(marker, region) %>%
  pull(region) %>%
  table() %>%
  prop.table()

all_markers_nearby_annotation %>%
  distinct(marker, region) %>%
  pull(region) %>%
  table() %>%
  prop.table()

# No difference in proportions




# Save the annotation data
save("all_sigmar_nearby_annotation", "all_markers_nearby_annotation", "egwas_sigmar", "eaa_gwas_spa_scores",
     "spa_outliers", "fst_outliers",  file = file.path(result_dir, "egwas_sigmar_analysis.RData"))























