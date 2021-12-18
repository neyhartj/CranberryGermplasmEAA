# Germplasm collection environmental association
#
# Figures and visualization
#
#

# Load packages
library(GenomicRanges)
library(raster)
library(tidyverse)
library(readxl)
library(broom)
library(neyhart)
library(patchwork)
library(cowplot)
library(scatterpie)
library(ggrepel)

# Load the startup script
source("startup.R")

# Set directories
cran_dir <- strsplit(proj_dir, "/")[[1]] %>%
  {.[seq_len(which(. == "CranberryLab"))]} %>%
  paste0(collapse = "/")

# Directory for environmental data
env_dir <- file.path(cran_dir, "EnvironmentalData")

# Labels for subfigures
subfigure_labels <- LETTERS

# Base font size
base_font_size <- 8

# Settings for figures
fig_dpi <- 300
fig_device <- "tiff"

# Create a new ggsave function with these defaults
ggsave2 <- function(filename, plot, device = fig_device, path = fig_dir,
                    height, width, units = c("in", "cm", "mm"), dpi = fig_dpi) {

  # Edit the filename to remove an extension
  if (str_detect(filename, "\\.")) {
    filename1 <- paste0(head(str_split(string = filename, pattern = "\\.")[[1]], -1), collapse = "_")
  } else {
    filename1 <- filename
  }
  filename2 <- paste0(filename1, ".", fig_device)

  # Use ggsave
  ggsave(filename = filename2, plot = plot, device = device, path = path, width = width,
         height = height, units = units, dpi = dpi)

}

# A function that replaces environmental variable short names with the full names
# given a width parameter
f_rename_variables <- function(x, width = 30, abbreviate = TRUE) {
  if (x %in% c("SPA", "F[ST]")){
    return(x)
  } else {
    full_names <- subset(eaa_environmental_vars, variable %in% x, full_name, drop = TRUE)
    if (abbreviate) {
      full_names <- full_names %>% str_replace_all(., "Temperature", "Temp.") %>%
        str_replace_all("Precipitation", "Precip.") %>%
        str_replace_all(., "Max", "Max.") %>%
        str_replace_all("Min", "Min.")
    }
    str_wrap(string = full_names, width = width)
  }
}



# Read in data ------------------------------------------------------------

# Read in the base data
load(file.path(data_dir, "population_metadata_and_genotypes.RData"))
load(file.path(data_dir, "germplasm_origin_bioclim_data.RData"))
load(file.path(result_dir, "population_genetics_stats.RData"))
# Load the eGWAS results
load(file.path(result_dir, "eaa_gwas_results.RData"))
load(file.path(result_dir, "egwas_sigmar_analysis.RData"))


# Load the annotation data
gff_file <- file.path(cran_dir, "Genotyping/ReferenceGenomes/Vaccinium_macrocarpon_BenLear_v2_annotations.gff")
# read in the GFF file
cranberry_gff <- ape::read.gff(file = gff_file)
cran_gene_ranges <- rtracklayer::import.gff3(con = gff_file)


# Calculate min/max chrom lengths
chromlengths <- cranberry_gff %>%
  filter(type == "chromosome", seqid != "Unknown") %>%
  select(chrom = seqid, start, end) %>%
  mutate(chrom = str_pad(chrom, 2, pad = "0")) %>%
  arrange(chrom)


# Prepare data ------------------------------------------------------------


# Select the location abbreviation and variables
eaa_environmental_data1 <- eaa_environmental_data %>%
  select(location_abbr, all_of(eaa_environmental_vars$variable))


# Prepare the bioclim data as if it were trait data
germplasm_bioclim_data <- pop_metadata %>%
  filter(category == "Wild", !is.na(latitude)) %>%
  select(individual, latitude, longitude) %>%
  left_join(., eaa_environmental_data1) %>%
  as.data.frame() %>%
  filter(individual %in% colnames(K_wild))

geno_mat_wild1 <- geno_mat_wild[germplasm_bioclim_data$individual,]

geno_hmp_wild1 <- geno_hmp_wild %>%
  select(marker:alleles, row.names(geno_mat_wild1))

# Tranpose the genotype hmp object
geno_hmp_wild2 <- geno_hmp_wild1 %>%
  separate(col = alleles, c("ref_allele", "alt_allele"), sep = "/") %>%
  mutate(marker1 = paste0(marker, "_", ref_allele, "_", alt_allele)) %>%
  filter(marker %in% unique(c(egwas_sigmar$marker, fst_outliers$marker, spa_outliers$marker))) %>%
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

# Data.frame for native selections and cultivars
geno_hmp_native_cultivars_wild <- geno_mat_all %>%
  {. + 1} %>% # Convert to ref allele count
  t() %>%
  as.data.frame() %>%
  rownames_to_column("marker") %>%
  inner_join(snp_info, .) %>%
  separate(alleles, c("ref_allele", "alt_allele"), sep = "/")


# Frequency of reference alleles for all markers in all population classes
all_marker_freq <- geno_hmp_native_cultivars_wild %>%
  gather(individual, ref_allele_count, -marker:-alt_allele) %>%
  left_join(., select(pop_metadata, individual, category)) %>%
  group_by(marker, ref_allele, category) %>%
  summarize(ref_allele_freq = mean(ref_allele_count, na.rm = TRUE) / 2, .groups = "drop")


# Marker metadata
marker_metadata <- geno_hmp_wild2 %>%
  distinct_at(vars(marker, chrom, pos, contains("allele")))


# Allele colors
allele_colors <- setNames(neyhart_palette("umn2")[3:4], c("ref_allele", "alt_allele"))


## Population sizes
# Summarize number of genotyped by site
wild_germ_meta_origin_summ <- pop_metadata %>%
  filter(category == "Wild") %>%
  rename(lat = latitude, long = longitude) %>%
  group_by(long, lat) %>%
  summarize(nSamples = n(), .groups = "drop")



## Plot sampling location

# Create a box for x and y axes
xlim <- range(pretty(germplasm_bioclim_data$long))
ylim <- range(pretty(germplasm_bioclim_data$lat))

# Map
g_map <- ggplot(data = north_america_mapdata, aes(x = long, y = lat)) +
  geom_polygon(fill = "white") +
  geom_polygon(data = subset(north_america_mapdata, area == "canada"),
               aes(group = group), fill = "grey95", color = "grey50", lwd = 0.5) + # Add canada
  geom_polygon(data = subset(north_america_mapdata, area == "usa_state"),
               aes(group = group), fill = "grey95", color = "grey50", lwd = 0.5) +
  coord_map(projection = "bonne", lat0 = mean(ylim), xlim = xlim, ylim = ylim) +
  theme_void(base_size = 14) +
  theme(legend.position = "top", legend.box = "vertical", legend.box.just = "left")








# Create lists of many plots ----------------------------------------------


### Create manhattan plots for each association ###

# Subset GWAS results for the best model
gwas_out_best_model <- gwas_out %>%
  # inner_join(., subset(p_lambda, best_model == "*", c(model, variable), drop = TRUE)
  filter(model == "model2")

# Manhattan plots
chrom_colors <- setNames(neyhart_palette("umn1")[3:4], c("odd", "even"))

eaa_gwas_manhattan_plot_df <- gwas_out %>%
  distinct(variable) %>%
  mutate(best_model = as.character(NA), plot = list(NULL))

# Tibble to list the significant markers for each variable
eaa_gwas_best_model_signif_marker <- eaa_gwas_manhattan_plot_df %>%
  mutate(sigmar = list(NULL)) %>%
  select(-plot)


# Iterate and plot - save each plot
for (i in seq_len(nrow(eaa_gwas_manhattan_plot_df))) {
  vari <- eaa_gwas_manhattan_plot_df$variable[i]

  scores1 <- gwas_out_best_model %>%
    filter(variable == vari) %>%
    mutate(chrom_color = ifelse(parse_number(chrom) %% 2 == 1, chrom_colors["odd"], chrom_colors["even"]))

  # Add the best model to the df
  selected_model <- unique(scores1$model)
  nPCs <- subset(egwas_models_df, model == selected_model, PCs, drop = TRUE)
  nPCs <- if (nPCs > 0) "Q" else character(0)
  selected_model_exp <- paste0(c("K", nPCs, subset(egwas_models_df, model == selected_model, covariates, drop = TRUE)[[1]]), collapse = " + ")
  eaa_gwas_manhattan_plot_df$best_model[i] <- selected_model

  # Full name of the vari
  vari_fullname <- subset(eaa_environmental_vars, variable == vari, full_name, drop = TRUE)

  # Get FDR threholds
  fdr_0.20 <- sommer:::fdr(p = scores1$p_value, fdr.level = 0.20)$fdr.10
  fdr_0.05 <- sommer:::fdr(p = scores1$p_value, fdr.level = 0.05)$fdr.10

  fdr_df <- tibble(annotation = paste0("FDR: ", c(0.05, 0.20)), fdr_p = c(fdr_0.05, fdr_0.20)) %>%
    mutate(fdr_score = -log10(fdr_p), annotation = fct_inorder(annotation))

  # Plot
  g_man <- scores1 %>%
    ggplot(aes(x = pos, y = score)) +
    geom_hline(data = fdr_df, aes(linetype = annotation, yintercept = fdr_score), color = "grey75") +
    geom_point(color = scores1$chrom_color, size = 1) +
    scale_x_continuous(breaks = median, name = NULL) +
    scale_y_continuous(name = expression(-log[10](p-value)), breaks = pretty, expand = expansion(mult = c(0, 0.35))) +
    scale_linetype_discrete(name = NULL) +
    facet_grid(~ chrom, scales = "free_x", space = "free_x", switch = "x") +
    labs(subtitle = vari_fullname,
         caption = paste0(str_to_title(selected_model), ": ", selected_model_exp)) +
    theme_genetics(base_size = 16) +
    theme(panel.spacing.x = unit(0, "lines"), strip.placement = "outside",
          axis.text.x = element_blank(), axis.ticks.length.x = unit(0.25, "lines"), axis.line.x = element_blank(),
          legend.position = c(0.99, 0.99), legend.justification = c(0.99, 0.99),
          legend.background = element_rect(fill = alpha("white", 0)),
          strip.background.x = element_blank())


  eaa_gwas_manhattan_plot_df$plot[[i]] <- g_man

  # Add significant markers to the df
  eaa_gwas_best_model_signif_marker$best_model[i] <- selected_model
  eaa_gwas_best_model_signif_marker$sigmar[[i]] <- scores1 %>%
    filter(p_value <= fdr_0.20) %>%
    select(marker, p_value)


  # filename <- paste0("egwas_manhattan_", vari, ".jpg")
  # ggsave(filename = filename, plot = g_man, path = fig_dir, width = 8, height = 4, dpi = 300)


}


# Sort the manhattan plots
eaa_gwas_manhattan_plot_df <- eaa_gwas_manhattan_plot_df %>%
  mutate(variable = factor(variable, levels = eaa_environmental_vars$variable)) %>%
  arrange(variable) %>%
  mutate(variable = as.character(variable))





# Iterate over markers and create two plots:
# 1. alleles versus bioclim variation
# 2. alleles across geography

# Combine all significant markers and FST/SPA outliers
all_sigmar <- egwas_sigmar %>%
  bind_rows(., select(mutate(fst_outliers, variable = "F[ST]"), marker, chrom, pos, score = Fst, variable)) %>%
  bind_rows(., select(mutate(spa_outliers, variable = "SPA"), marker, chrom, pos, score = spa_score, variable))

#
sigmar_summary <- all_sigmar %>%
  bind_rows(., select(filter(gwas_out_best_model, marker == "S07_30582060", variable == "bio1"), variable, best_model = model, marker, p_value)) %>%
  group_by(variable, marker) %>%
  do({
    row <- .
    marker_i <- row$marker
    var_i <- row$variable
    var_i_fullname <- subset(eaa_environmental_vars, variable == var_i, full_name, drop = TRUE)


    # Biovar selector depending on variable
    if (var_i %in% eaa_environmental_vars$variable) {
      # Get the mlmm output for this variable
      biovar1 <- select(germplasm_bioclim_data, individual, biovar1 = all_of(var_i))
      location_biovar <- select(eaa_environmental_data, location_abbr, biovar = all_of(var_i))

    } else {
      biovar1 <- select(germplasm_bioclim_data, individual) %>% mutate(biovar1 = as.numeric(NA))
      location_biovar <- select(eaa_environmental_data, location_abbr)

    }


    # Subset the bioclim data for this variable
    bioclim_data_i <- germplasm_bioclim_data %>%
      select(individual, lat = latitude, long = longitude, location_abbr) %>%
      left_join(., location_biovar, by = "location_abbr") %>%
      # Add the SNP alleles
      left_join(., select(geno_hmp_wild1t, individual, genotype = all_of(marker_i)), by = "individual") %>%
      group_by(lat, long, genotype) %>%
      mutate(genotypes_by_location = n()) %>%
      group_by(genotype) %>%
      mutate(genotypes_overall = n()) %>%
      ungroup() %>%
      mutate_at(vars(genotypes_by_location, genotypes_overall), ~paste0(genotype, "\n(n=", ., ")")) %>%
      left_join(., biovar1, by = "individual")

    # Subset reference and alternate alleles for the marker
    ref_allele <- subset(marker_metadata, marker == marker_i, ref_allele, drop = TRUE)
    alt_allele <- subset(marker_metadata, marker == marker_i, alt_allele, drop = TRUE)


    # Recreate but with ref allele counts
    bioclim_data_i1 <- germplasm_bioclim_data %>%
      select(individual, lat = latitude, long = longitude, location_abbr) %>%
      left_join(., location_biovar, by = "location_abbr") %>%
      # Add the SNP alleles
      left_join(., select(geno_hmp_wild1t2, individual, genotype = all_of(marker_i)), by = "individual") %>%
      group_by(lat, long, location_abbr) %>%
      summarize(freq = mean(genotype) / 2, .groups = "drop") %>%
      mutate(allele = ref_allele)

    bioclim_data_i1 <- bioclim_data_i1 %>%
      bind_rows(., mutate(bioclim_data_i1, allele = alt_allele,
                          freq = 1 - freq))

    # Spread out alleles to columns
    bioclim_data_i2 <- bioclim_data_i1 %>%
      spread(allele, freq) %>%
      # add population size
      left_join(., wild_germ_meta_origin_summ, by = c("lat", "long"))



    # Create genotype colors
    genotype_colors <- colorRampPalette(colors = allele_colors)(3) %>%
      setNames(., c(paste0(ref_allele, ref_allele), paste0(ref_allele, alt_allele), paste0(alt_allele, alt_allele)))


    # Plot marker allele versus bioclim
    g1 <- ggplot(data = bioclim_data_i, aes(x = genotypes_overall, y = biovar1, group = genotype)) +
      # geom_boxplot(width = 0.5, alpha = 0.5, outlier.shape = NA) +
      geom_violin(aes(fill = genotype), alpha = 0.5) +
      geom_jitter(width = 0.1, size = 0.5) +
      xlab(marker_i) +
      ylab(str_wrap(var_i_fullname, width = 30)) +
      scale_fill_manual(values = genotype_colors, guide = FALSE) +
      neyhart::theme_genetics()

    colors <- structure(RColorBrewer::brewer.pal(n = 11, name = "RdBu"), class = "palette")[c(3, 6, 9)]
    colors <- RColorBrewer::brewer.pal(n = 9, name = "Blues")[c(2, 5, 8)]

    # Vector of radii
    radii <- (bioclim_data_i2$nSamples^0.33) / 5
    # Function to convert
    f_revert_radii <- function(x) (x * 5)^(1/0.33)


    # Rename allele colors
    allele_colors1 <- allele_colors
    names(allele_colors1) <- c(ref_allele, alt_allele)

    # Plot frequency of alleles in different geographies
    # Plot marker alleles across geographies
    suppressMessages({
      g2 <- g_map +
        # Scatterpie
        geom_scatterpie(data = bioclim_data_i2, aes(x = long, y = lat, r = (nSamples^0.33) / 5),
                        cols = c(ref_allele, alt_allele), alpha = 0.65) +
        coord_fixed(xlim = xlim, ylim = ylim, ratio = 1) +
        # geom_scatterpie_legend(radius = radii, x = max(xlim) - 2, y = min(ylim) + 2, labeller = f_revert_radii, n = 2)
        scale_fill_manual(values = allele_colors1, name = paste0(marker_i, "\nmarker allele"),
                          guide = guide_legend(title.position = "top")) +
        theme(legend.position = c(0.8, 0.2), legend.margin = margin(3,3,3,3), legend.direction = "horizontal")
    })

    # Alternative map with individual points and intermediate colors for hets
    g2_alt <- g_map +
      geom_jitter(data = bioclim_data_i, aes(color = genotype), width = 0.25, height = 0.25, size = 1) +
      # scale_color_manual(values = genotype_colors, name = "Marker\ngenotype") +
      # theme(legend.position = c(0, 0), legend.justification = c(0, 0),
      #       legend.box.background = element_rect(fill = "white", colour = NA), legend.box.margin = margin(2, 5, 2, 5))
      scale_color_manual(values = genotype_colors, name = "Marker genotype",
                         guide = guide_legend(override.aes = list(size = 2))) +
      theme(legend.position = "bottom", legend.justification = "left")

    # Return a tibble
    tibble(plot1 = list(g1), plot2 = list(g2), plot2_alt = list(g2_alt))


  }) %>% ungroup() %>%
  left_join(., select(eaa_environmental_vars, -class)) %>%
  mutate(variable = factor(variable, levels = eaa_environmental_vars$variable)) %>%
  arrange(variable)



# # Create a bunch of plots
# pb <- progress::progress_bar$new(total = nrow(sigmar_summary))
#
# for (i in seq_len(nrow(sigmar_summary))) {
#
#   # Determine the left and right plots
#   left_plot <- sigmar_summary$plot2[[i]] +
#     theme(legend.position = "none")
#   right_plot <- sigmar_summary$plot1[[i]]
#
#   combined_plot1 <- plot_grid(left_plot, right_plot, rel_widths = c(1, 0.4))
#
#   # Save
#   filename <- paste0("egwas_summary_", sigmar_summary$variable[i], "-marker-", sigmar_summary$marker[i], "-",
#                      sigmar_summary$class[i], ".jpg")
#   ggsave(filename = filename, plot = combined_plot1, path = fig_dir, width = 8, height = 2.5, dpi = 300)
#
#   pb$tick()
#
# }





# Figure 1: geographic origin, population structure, ld decay -------------

# Plot origin of wild germplasm

# Create a box for x and y axes
xlim <- range(pretty(wild_germ_meta1$long))
ylim <- range(pretty(wild_germ_meta1$lat))

# Breaks for population size
pop_size_breaks <- pretty(range(wild_germ_meta_origin_summ$nSamples))

# Map
g_map <- ggplot(data = north_america_mapdata, aes(x = long, y = lat)) +
  geom_polygon(fill = "white") +
  geom_polygon(data = subset(north_america_mapdata, area == "canada"),
               aes(group = group), fill = "grey95", color = "grey50", lwd = 0.5) + # Add canada
  geom_polygon(data = subset(north_america_mapdata, area == "usa_state"),
               aes(group = group), fill = "grey95", color = "grey50", lwd = 0.5) +
  geom_point(data = wild_germ_meta_origin_summ, aes(size = nSamples)) +
  coord_map(projection = "bonne", lat0 = mean(ylim), xlim = xlim, ylim = ylim) +
  scale_size_continuous(name = "Population size", breaks = pop_size_breaks, labels = pop_size_breaks,
                        limits = range(pop_size_breaks), range = c(1, 8)) +
  # scale_size_manual() +
  labs(subtitle = "Wild cranberry collection sites") +
  theme_void(base_size = 14) +
  theme(legend.position = c(0.90, 0.05), legend.justification = c(1,0), legend.box = "vertical", legend.box.just = "left",
        panel.border = element_rect(fill = NA))


# Save the figure
ggsave(filename = "cranberry_wild_germplasm_origins.jpg", plot = g_map, path = fig_dir,
       width = 8, height = 4, dpi = 1000)

# Remove the subtitle
g_map1 <- g_map
g_map1$labels$subtitle <- NULL


# Save the figure
ggsave2(filename = "figure1_cranberry_wild_germplasm_origins1", plot = g_map1, path = fig_dir,
        device = fig_device, width = 8, height = 4, dpi = fig_dpi)





# Figure 2: Population structure and LD decay --------------------------------------


# PCA of wild individuals for determining PCs
K_pca <- prcomp(x = K_wild)
# PCA of all accessions
K_pca_all <- prcomp(x = K_all)

# Gather the eigenvector
pca_tidy <- broom::tidy(K_pca) %>%
  rename(individual = row) %>%
  mutate(PC = paste0("PC", PC)) %>%
  spread(PC, value) %>%
  left_join(., select(pop_metadata, individual, category, location_of_origin)) %>%
  left_join(., select(eaa_environmental_data, location_of_origin, location_abbr, state))

# Get the eigenvalues
eigenvals <- K_pca$sdev^2

# Calculate proportion of explained variance
pc_varexp <- eigenvals / sum(eigenvals)
names(pc_varexp) <- paste0("PC", seq_along(eigenvals))

# Cumulative sum
plot(cumsum(pc_varexp))
head(pc_varexp)
sum(pc_varexp[1:3])

# Calculate the median position of each color and plot a label for that
group_label <- pca_tidy %>%
  group_by(state) %>%
  summarize_at(vars(contains("PC")), median)

# add Origin and plot
pc_plot <- pca_tidy %>%
  ggplot(aes(x = PC1, y = PC2, color = state)) +
  geom_point(aes(alpha = location_abbr, linetype = individual)) +
  ggrepel::geom_label_repel(data = group_label, aes(label = state), fill = NA, key_glyph = "point", size = 2,
                            max.overlaps = 15, point.padding = unit(10, "line"), box.padding = unit(0.2, "line")) +
  scale_color_discrete(guide = FALSE) +
  scale_x_continuous(name = paste0("PC1 (", round(pc_varexp["PC1"] * 100, 3), "%)"), breaks = pretty,
                     labels = format_numbers) +
  scale_y_continuous(name = paste0("PC2 (", round(pc_varexp["PC2"] * 100, 3), "%)"), breaks = pretty,
                     labels = format_numbers) +
  scale_alpha_discrete(guide = FALSE) +
  theme_genetics() +
  theme(legend.position = c(0.99, 0.99), legend.justification = c(0.99, 0.99))

# Save this
ggsave(filename = "snp_pca_wild_germplasm.jpg", plot = pc_plot, path = fig_dir, width = 5, height = 4, dpi = 1000)


## Plot PCA of all entries

# Gather the eigenvector
pca_tidy <- broom::tidy(K_pca_all) %>%
  rename(individual = row) %>%
  mutate(PC = paste0("PC", PC)) %>%
  spread(PC, value) %>%
  left_join(., select(pop_metadata, individual, category, location_of_origin)) %>%
  left_join(., select(eaa_environmental_data, location_of_origin, location_abbr, state)) %>%
  # Assign no location to NA
  mutate(state = ifelse(is.na(state), "(none)", state))

# Get the eigenvalues
eigenvals <- K_pca_all$sdev^2

# Calculate proportion of explained variance
pc_varexp <- eigenvals / sum(eigenvals)
names(pc_varexp) <- paste0("PC", seq_along(eigenvals))

# Cumulative sum
plot(cumsum(pc_varexp))
head(pc_varexp)
sum(pc_varexp[1:3])

# Calculate the median position of each color and plot a label for that
group_label <- pca_tidy %>%
  group_by(state) %>%
  summarize_at(vars(contains("PC")), median)

# add Origin and plot
pc_plot_all <- pca_tidy %>%
  ggplot(aes(x = PC1, y = PC2, color = state, shape = category)) +
  geom_point() +
  scale_color_manual(values = c("black", RColorBrewer::brewer.pal(n = 12, name = "Set3")[-2])) +
  scale_x_continuous(name = paste0("PC1 (", round(pc_varexp["PC1"] * 100, 3), "%)"), breaks = pretty,
                     labels = format_numbers) +
  scale_y_continuous(name = paste0("PC2 (", round(pc_varexp["PC2"] * 100, 3), "%)"), breaks = pretty,
                     labels = format_numbers) +
  scale_alpha_discrete(guide = FALSE) +
  theme_genetics() +
  theme(legend.position = "right", legend.justification = c(0.99, 0.99))

# Save this
ggsave(filename = "snp_pca_all_germplasm.jpg", plot = pc_plot_all, path = fig_dir, width = 6.5, height = 5, dpi = 1000)




###

# Modify the fitted LD values
ld_wild_decay1 <- unnest(ld_wild_decay, predictions)
ld_wild_df1 <- ld_wild_df %>%
  filter(bp_dist <= max(ld_wild_decay1$d))

## Plot LD decay

# Plot LD decay without filtering distance

g_ld_decay_all_chrom1 <- ld_wild_df1 %>%
  ggplot(aes(x = bp_dist / 1000, y = r2)) +
  geom_point(size = 0.3, color = "gray") +
  geom_line(data = ld_wild_decay1, aes(x = d / 1000, y = yhat), lwd = 1) +
  scale_x_continuous(name = "Distance (kbp)", breaks = pretty) +
  scale_y_continuous(name = expression("Linkage disequilibrium ("*r^2*")")) +
  facet_wrap(~ LG, labeller = labeller(LG = function(x) paste0("Chrom. ", x))) +
  theme_genetics()


# Plot LD decay with filtering distance

# Modify the fitted LD values
ld_wild_decay2 <- unnest(ld_wild_decay, predictions) %>%
  filter(d <= 2e5)
ld_wild_df2 <- ld_wild_df %>%
  # rename(d = bp_dist) %>%
  filter(bp_dist <= 2e5)


g_ld_decay_all_chrom2 <- ld_wild_df2 %>%
  ggplot(aes(x = bp_dist / 1000, y = r2)) +
  geom_point(size = 0.3, color = "gray") +
  geom_line(data = ld_wild_decay2, aes(x = d / 1000, y = yhat), lwd = 1) +
  scale_x_continuous(name = "Distance (kbp)", breaks = pretty) +
  scale_y_continuous(name = expression("Linkage disequilibrium ("*r^2*")")) +
  facet_wrap(~ LG, labeller = labeller(LG = function(x) paste0("Chrom. ", x))) +
  theme_genetics()


## Plot chromosome 1 low distance versus full

g_ld_decay_chrom1_long <- ld_wild_df1 %>%
  filter(LG == "01") %>%
  ggplot(aes(x = bp_dist / 1000, y = r2)) +
  geom_point(size = 0.3, color = "gray") +
  geom_line(data = subset(ld_wild_decay1, LG == "01"), aes(x = d / 1000, y = yhat), lwd = 1) +
  scale_x_continuous(name = "Chromosome 1 marker distance (kbp)", breaks = pretty, expand = expansion(mult = c(0.01, 0.05))) +
  scale_y_continuous(name = expression("Linkage disequilibrium ("*r^2*")"), expand = expansion(mult = c(0.01, 0.05))) +
  theme_genetics()

g_ld_decay_chrom1_short <- ld_wild_df2 %>%
  filter(LG == "01") %>%
  ggplot(aes(x = bp_dist / 1000, y = r2)) +
  geom_point(size = 0.3, color = "gray") +
  geom_line(data = subset(ld_wild_decay2, LG == "01"), aes(x = d / 1000, y = yhat), lwd = 1) +
  scale_x_continuous(name = "Chromosome 1 marker distance (kbp)", breaks = pretty, expand = expansion(mult = c(0.01, 0.01))) +
  scale_y_continuous(name = expression("Linkage disequilibrium ("*r^2*")"), expand = expansion(mult = c(0.01, 0.05))) +
  theme_genetics()



# Determine the left and right plots
left_plot <- g_ld_decay_chrom1_short
right_plot <- g_ld_decay_chrom1_long +
  theme(axis.title = element_blank(), plot.background = element_rect(fill = alpha("white", 1), color = "black"))

## Combine plots
layout <- c(
  area(t = 1, l = 1, b = 10, r = 8),
  area(t = 1, l = 6, b = 4, r = 8)
)

combined_ld_plot <- plot_grid(left_plot) +
  plot_grid(right_plot)+
  plot_layout(design = layout)



# Combine pop structure with LD
g_popstr_ld <- plot_grid(pc_plot, combined_ld_plot, ncol = 1, labels = subfigure_labels[1:2],
                         label_size = 9, align = "h", axis = "tblr")

# Save
ggsave2(filename = "figure2_popstr_chrom1_ld", plot = g_popstr_ld, width = 3.5, height = 6, dpi = 500)



#
#
#
#
# # What is the average pairwise genetic distance between populations?
# pairwise_dist_by_loc <- pairwise_dist_data %>%
#   filter(location1 != location2) %>%
#   mutate(locs_merge = map2_chr(location1, location2, ~paste0(sort(c(.x, .y)), collapse = "_"))) %>%
#   group_by(locs_merge) %>%
#   summarize(mean_distance = mean(distance), .groups = "drop") %>%
#   separate(locs_merge, c("location1", "location2"), sep = "_") %>%
#   spread(location1, mean_distance) %>%
#   as.data.frame() %>%
#   column_to_rownames("location2") %>%
#   as.dist()
#
# # Calculate average relationship between individuals in different locations
# pairwise_fXY_by_loc <- K_wild %>%
#   {. / 2} %>% # Off-diagonal elements should be 2 x f_XY
#   as.data.frame() %>%
#   rownames_to_column("individual1") %>%
#   gather(individual2, f_XY, -individual1) %>%
#   left_join(., distinct(germplasm_bioclim_data, individual, location_abbr), by = c("individual1" = "individual")) %>%
#   rename(location1 = location_abbr) %>%
#   left_join(., distinct(germplasm_bioclim_data, individual, location_abbr), by = c("individual2" = "individual")) %>%
#   rename(location2 = location_abbr) %>%
#   filter(location1 != location2) %>%
#   mutate(locs_merge = map2_chr(location1, location2, ~paste0(sort(c(.x, .y)), collapse = "_"))) %>%
#   group_by(locs_merge) %>%
#   summarize(f_XY = mean(f_XY), .groups = "drop") %>%
#   separate(locs_merge, c("location1", "location2"), sep = "_") %>%
#   spread(location1, f_XY) %>%
#   as.data.frame() %>%
#   column_to_rownames("location2") %>%
#   as.dist()




# Figure 3: Summary of SPA, Fst, and EGWAS --------------------------------

# Combine results for all three analyses for a single image
all_stats_to_plot <- egwas_sigmar %>%
  select(marker, chrom, pos, variable) %>%
  left_join(., eaa_environmental_vars) %>%
  bind_rows(mutate(select(fst_outliers, -Fst, -Fis), variable = "F[ST]", full_name = variable, class = variable)) %>%
  bind_rows(mutate(select(spa_outliers, chrom, marker, pos), variable = "SPA", full_name = variable, class = variable)) %>%
  group_by(marker, class) %>%
  mutate(nVar = n()) %>%
  ungroup() %>%
  filter(class != "principal_component") %>%
  mutate(nVar = ifelse(class %in% c("F[ST]", "SPA"), max(nVar), nVar),
         # rename class
         class = case_when(
           class == "soil" ~ "EAA~(Soil)",
           class == "temperature" ~ "EAA~(Temp.)",
           class == "precipitation" ~ "EAA~(Prec.)",
           class == "geography" ~ "EAA~(Geo.)",
           TRUE ~ as.character(class)),
         class = fct_inorder(class) %>% fct_rev())

# Unique
all_stats_to_plot_unique <- all_stats_to_plot %>%
  distinct(marker, chrom, pos, class, nVar)


# First a base plot so all chromosome lengths are accounted
g_all_stats <- chromlengths %>%
  gather(mark, pos, -chrom) %>%
  ggplot(aes(x = pos / 1e6)) +
  geom_point(data = all_stats_to_plot, aes(y = class, alpha = as.factor(nVar)), pch = 16, color = "red") +
  facet_wrap(~ chrom, ncol = 4, labeller = labeller(chrom = function(x) paste0("Chrom. ", x))) +
  scale_x_continuous(name = "Position (Mbp)", breaks = pretty) +
  scale_y_discrete(name = NULL, labels = function(x) parse(text = x)) +
  scale_alpha_manual(guide = FALSE, values = seq(0.5, 1, length.out = 6)) +
  theme_genetics() +
  theme(panel.border = element_rect(fill = alpha("white", 0)), strip.background = element_blank(),
        panel.grid.major.y = element_line(color = "gray85"), panel.grid.major.x = element_line(color = "gray85"),
        panel.grid.minor.x = element_line(color = "gray85"))

# Save
ggsave2(filename = "figure3_genomic_statistics", plot = g_all_stats, height = 3, width = 8, dpi = 500)



## Redesign using chromosome shapes

chromlengths1 <- mutate(chromlengths, start = 0)
chromlengths2 <- gather(chromlengths1, position, bp, -chrom) %>%
  mutate(pos_ann = ifelse(bp != 0, paste0(round(bp / 1e6), "\nMb"), bp))

# Assign class labels to SNPs within a certain distance
all_stats_to_plot_unique_labels <- all_stats_to_plot_unique %>%
  arrange(chrom, pos) %>%
  group_by(chrom) %>%
  do({
    df <- .
    df1 <- mutate(df, pos_diff = diff(c(0, pos)))
    which_label <- union(which(df1$pos_diff < 1e6), which(df1$pos_diff < 1e6)-1)
    df1$label = as.character(NA)
    df1$label[which_label] <- as.character(df1$class[which_label])
    df1
  }) %>% ungroup()

# Segment width adjuster
segWidth <- 200000

# First create a plot of blank chromosomes
g_chromlengths <- chromlengths1 %>%
  ggplot(aes(y = start, yend = end, x = as.factor(0), xend = as.factor(0))) +
  geom_segment(lwd = 5, lineend = "round", color = "grey85") +
  # Add start/end text
  geom_text(data = chromlengths2, aes(x = as.factor(0), y = bp, label = pos_ann), inherit.aes = FALSE,
            hjust = 1, nudge_x = -0.2, size = 3) +
  geom_segment(data = chromlengths2, aes(y = bp - segWidth, yend = bp + segWidth), lwd = 5, color = "black") +
  # Add annotations
  geom_segment(data = subset(all_stats_to_plot_unique),
               aes(y = pos - segWidth, yend = pos + segWidth, color  = class, alpha = as.factor(nVar)), lwd = 5) +
  geom_text_repel(data = all_stats_to_plot_unique_labels, aes(x = as.factor(0), y = pos, label = label, color = class),
                  inherit.aes = FALSE, direction = "y", size = 2, hjust = 0, parse = TRUE,
                  nudge_x = 4, force_pull = 0, force = 0.25, min.segment.length = 0) +
  facet_wrap(~ chrom, strip.position = "bottom", nrow = 2) +
  scale_alpha_manual(guide = FALSE, values = seq(0.5, 1, length.out = 6)) +
  scale_color_manual(values = c(neyhart_palette("umn2")[3], neyhart_palette("fall")), name = NULL,
                     guide = guide_legend(nrow = 1), labels = function(x) parse(text = x)) +
  theme_nothing(font_size = base_font_size) +
  theme(strip.text.x = element_text(size = base_font_size), panel.spacing = unit(0, "line"),
        legend.position = "top")

ggsave2(filename = "figure3_genomic_statistics_chroms", plot = g_chromlengths,
        width = 8, height = 5.5, dpi = 500)




# Figure 4: cold temperature association example -----------------------------


## First show the locus on Chromosome 4 associated with multiple worldclim temperature
## variables

# loci_highlight <- c("S04_24963382", "S07_30582060")
loci_highlight <- c("S04_24963382")

egwas_sigmar %>%
  filter(marker %in% loci_highlight) %>%
  # add environmental variable full names
  left_join(., select(eaa_environmental_vars, -class))

locus1_summary <- egwas_sigmar %>%
  filter(marker %in% loci_highlight, str_detect(variable, "bio")) %>%
  # add environmental variable full names
  left_join(., select(eaa_environmental_vars, -class))

# Create a highlight window for this snp
locus1_highlight_rect <- snp_info %>%
  filter(marker %in% loci_highlight) %>%
  mutate(left = pos - 1000, right = pos + 1000)


## Plot Part A - combined manhattan plots

# Get the manhattan plot data for all of the variable associated with the loci
manhattan_data <- eaa_gwas_manhattan_plot_df %>%
  filter(variable %in% locus1_summary$variable) %>%
  mutate(plot_data = map(plot, "data"), variable = factor(variable, levels = eaa_environmental_vars$variable)) %>%
  select(variable, plot_data) %>%
  unnest(plot_data, names_repair = tidyr_legacy)

# Create a DF for adding fdr thresholds
fdr_df <- manhattan_data %>%
  group_by(variable) %>%
  nest() %>%
  # Calculate FDR levels for plotting
  mutate(fdr_0.20 = map_dbl(data, ~sommer:::fdr(p = .x$p_value, fdr.level = 0.20)$fdr.10),
       fdr_0.05 = map_dbl(data, ~sommer:::fdr(p = .x$p_value, fdr.level = 0.05)$fdr.10)) %>%
  mutate_at(vars(contains("fdr")), ~-log10(.)) %>%
  ungroup() %>%
  select(-data) %>%
  gather(fdr_level, fdr_score, -variable) %>%
  mutate(annotation = str_to_upper(str_replace_all(fdr_level, "_", " ")),
         annotation = fct_relevel(annotation, "FDR 0.05"))


# Plot
g_man_combined <- manhattan_data %>%
  ggplot(aes(x = pos, y = score)) +
  geom_hline(data = fdr_df, aes(linetype = annotation, yintercept = fdr_score), color = "grey75") +
  geom_point(color = manhattan_data$chrom_color, size = 0.5) +
  # Add highlighting line
  geom_vline(data = locus1_highlight_rect, aes(xintercept = pos), lwd = 2.5, color = "grey75", alpha = 0.2) +
  facet_grid(variable ~ chrom, scales = "free_x", space = "free_x", switch = "both",
             labeller = labeller(variable = function(x) f_rename_variables(x, width = 18))) +
  scale_x_continuous(breaks = median, name = NULL) +
  scale_y_continuous(name = expression(-log[10](p-value)), breaks = pretty, expand = expansion(mult = c(0, 0.35))) +
  scale_linetype_discrete(name = NULL) +
  theme_genetics(base_size = base_font_size) +
  theme(panel.spacing.x = unit(0, "lines"), strip.placement = "outside",
        axis.text.x = element_blank(), axis.ticks.length.x = unit(0.25, "lines"), axis.line.x = element_blank(),
        legend.position = c(0.99, 0.99), legend.justification = c(0.99, 0.99), legend.direction = "horizontal",
        legend.background = element_rect(fill = alpha("white", 0)),
        strip.background.x = element_blank())


## Plot Part B - Genomic range of the gene

loci_grange <- subset(snp_info, marker %in% loci_highlight) %>%
  mutate(start = pos - 17000, end = pos + 17000, chrom = parse_number(chrom)) %>%
  as(., "GRanges")

# Nearby overlaps
loci_overlaps <- subsetByOverlaps(cran_gene_ranges, loci_grange) %>%
  subset(type == "gene")

# For each gene, get exons
loci_overlaps_exons <- loci_overlaps %>%
  split(.$ID) %>%
  lapply(X = ., FUN = function(x) subsetByOverlaps(cran_gene_ranges, x)) %>%
  map(~subset(., type == "exon")) %>%
  imap_dfr(~as.data.frame(.x) %>% mutate(.id = .y)) %>%
  mutate(y = as.numeric(as.factor(strand)))

# Make introns
loci_overlaps_introns <- loci_overlaps_exons %>%
  split(.$.id) %>%
  map_df(~{
    mutate(head(.x, -1), start1 = end, end1 = c(start[-1], tail(.x$start, 1))) %>%
      split(.$ID) %>%
      map_df(~{
        mid <- floor(mean(c(.x$start1[1], .x$end1[1])))
        data.frame(ID = .x$ID, strand = .x$strand, start1 = c(.x$start1[1], mid), end1 = c(mid, .x$end1[1]),
                   y = c(.x$y[1], .x$y + 0.05), yend = c(.x$y + 0.05, .x$y[1]))
      })
  })



loci_overlaps_exons_range <- loci_overlaps_exons %>%
  group_by(.id, strand, y) %>%
  summarize(start = min(start), end = max(end), .groups = "drop")

loci_overlaps_df <- as.data.frame(loci_overlaps)

# Visualize
g_loci_overlaps <- loci_overlaps_exons %>%
  ggplot(aes(x = start, xend = end, y = y, yend = y, color = strand)) +
  # geom_segment(data = loci_overlaps_exons_range, lwd = 0.5) +
  # geom_segment(lwd = 4, arrow = arrow(length = unit(0.05, "line"), ends = ifelse(loci_overlaps_df$strand == "+", "first", "last"),
  #                                     angle = 2, type = "closed"), linejoin = "mitre") +
  geom_segment(lwd = 1.5, lineend = "square") +
  geom_segment(data = loci_overlaps_introns, aes(x = start1, xend = end1, yend = yend)) +
  geom_vline(data = locus1_highlight_rect, aes(xintercept = pos), lwd = 2.5, color = "grey75", alpha = 0.2) +
  scale_color_discrete(guide = FALSE) +
  # Adjust 'expand' argument to align the vertical lines from the manhattan plot and the gene model plot
  scale_x_continuous(name = "Chromosome 4 position", breaks = pretty) +
  scale_y_continuous(expand = expansion(mult = 0.5)) +
  theme_genetics(base_size = base_font_size) +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.line.y = element_blank())

# Offset amount
offset <- 0.10

# Polygons
loci_overlaps_polygons <- loci_overlaps_df %>%
  select(strand, start, end, Name) %>%
  mutate(strand_num = as.numeric(as.factor(strand))) %>%
  split(.$Name) %>%
  map_df(~{
    rect <- data.frame(strand = .x$strand,
               x = c(rep(.x$start, 2), rep(.x$end, 2)),
               y = c(c(.x$strand_num - offset, .x$strand_num + offset), rev(c(.x$strand_num - offset, .x$strand_num + offset))),
               Name = .x$Name,
               group = paste0(.x$Name, "_rect"))
    # Triangle
    xend <- ifelse(.x$strand == "+", .x$start - 1000, .x$end + 1000)
    xstart <- ifelse(.x$strand == "+", .x$start, .x$end)
    tri <- data.frame(strand = .x$strand,
                      x = c(xstart, xstart, xend),
                      y = c(.x$strand_num - offset, .x$strand_num + offset, .x$strand_num),
                      Name = .x$Name,
                      group = paste0(.x$Name, "_tri"))
    # Merge
    bind_rows(rect, tri)
  })

# DF for labels
loci_overlaps_labels <- loci_overlaps_polygons %>%
  group_by(Name) %>%
  summarize(x = (max(x) + min(x))/2, y = min(y) - (2*offset), strand = first(strand))

g_loci_overlaps <- loci_overlaps_polygons %>%
  ggplot(aes(x = x / 1e3, y = y)) +
  geom_polygon(aes(fill = strand, group = group)) +
  geom_vline(data = locus1_highlight_rect, aes(xintercept = pos / 1e3), lwd = 2.5, color = "grey75", alpha = 0.2) +
  geom_text(data = loci_overlaps_labels, mapping = aes(label = Name), size = 2, hjust = "inward") +
  # geom_text_repel(data = loci_overlaps_labels, mapping = aes(label = Name), size = 2, direction = "y") +
  scale_fill_discrete(guide = FALSE) +
  # Adjust 'expand' argument to align the vertical lines from the manhattan plot and the gene model plot
  scale_x_continuous(name = "Chromosome 4 position (kbp)", breaks = pretty) +
  scale_y_continuous(expand = expansion(mult = 0.1)) +
  theme_genetics(base_size = base_font_size) +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.line.y = element_blank())



## Part C - geographic distribution of the marker alleles

# Subset this from the list of plots above
allele_geo_list <- sigmar_summary %>%
  filter(marker %in% loci_highlight) %>%
  split(.$marker) %>%
  map("plot2") %>%
  map(1) %>%
  map(~.x + theme_nothing(font_size = base_font_size) + theme(legend.position = .x$theme$legend.position))


## Part D - Bioclim variables versus marker alleles

bioclim_allele_list <- sigmar_summary %>%
  filter(marker %in% loci_highlight, str_detect(variable, "bio")) %>%
  split(.$marker) %>%
  map(~{
    # Combine plot data
    plot_data_combined <- .x %>%
      mutate(data = map(plot1, "data")) %>%
      select(-contains("plot")) %>%
      unnest(data) %>%
      mutate(genotypes_overall = genotype)

    alleles <- str_split(subset(snp_info, marker == unique(.x$marker), alleles, drop = TRUE), "/")[[1]]
    ref_allele <- alleles[1]
    alt_allele <- alleles[2]

    # Create genotype colors
    genotype_colors <- colorRampPalette(colors = allele_colors)(3) %>%
      setNames(., c(paste0(ref_allele, ref_allele), paste0(ref_allele, alt_allele), paste0(alt_allele, alt_allele)))

    # Means by genotype class
    genotype_means <- aggregate(biovar ~ genotypes_overall + variable, plot_data_combined, mean)

    # Create a new plot and return
    plot_data_combined %>%
      ggplot(aes(x = 0, y = biovar, fill = genotypes_overall, group = genotypes_overall)) +
      geom_jitter(pch = 21, position = position_jitterdodge(jitter.width = 0.10, dodge.width = 0.75), size = 1) +
      geom_point(data = genotype_means, position = position_dodge(0.75), pch = "+", size = 5, color = neyhart_palette()[3]) +
      facet_wrap(~ variable, ncol = 3, scales = "free_y", dir = "v",
                 labeller = labeller(variable = f_rename_variables)) + # dir = "v" means fill by column
      scale_fill_manual(values = genotype_colors, guide = FALSE) +
      scale_x_continuous(name = NULL, labels = NULL, breaks = NULL) +
      scale_y_continuous(name = NULL, breaks = pretty) +
      theme_grey(base_size = base_font_size) +
      theme(strip.background = element_blank(), panel.grid.minor.x = element_blank(),
            panel.grid.major.x = element_blank())

  })


## Combine everything

p1 <- g_man_combined
p2 <- g_loci_overlaps
p3 <- allele_geo_list$S04_24963382
p4 <- bioclim_allele_list$S04_24963382

# For patchwork, you need to specify the nesting design, as demonstrated below
g_combined <- ( (p1 + p2 + plot_layout(ncol = 1, heights = c(1, 0.12))) |
    (p3 + p4 + plot_layout(ncol = 1, heights = c(1, 1))) ) +
  plot_layout(ncol = 2, widths = c(0.40, 1)) +
  plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(size = base_font_size + 2, face = "plain"))


# Save
ggsave2(filename = "figure4_chrom4_temperature_association", plot = g_combined,
        width = 8, height = 6, dpi = 500)



# What is the difference in variables between those carrying the minor/major allele at
# each marker?
bioclim_allele_list$S04_24963382$data %>%
  group_by(variable, genotype) %>%
  summarize(biovar_mean = mean(biovar), nInd = n(), .groups = "drop") %>%
  filter(map_lgl(str_split(genotype, ""), ~n_distinct(.) == 1)) %>%
  mutate(genotype_class = ifelse(nInd == min(nInd), "minor", "major")) %>%
  select(-genotype, -nInd) %>%
  spread(genotype_class, biovar_mean) %>%
  mutate(diff = minor - major)

bioclim_allele_list$S07_30582060$data %>%
  group_by(variable, genotype) %>%
  summarize(biovar_mean = mean(biovar), nInd = n(), .groups = "drop") %>%
  filter(map_lgl(str_split(genotype, ""), ~n_distinct(.) == 1)) %>%
  mutate(genotype_class = ifelse(nInd == min(nInd), "minor", "major")) %>%
  select(-genotype, -nInd) %>%
  spread(genotype_class, biovar_mean) %>%
  mutate(diff = minor - major)


# Frequency of alleles in all populations
subset(all_marker_freq, marker %in% loci_highlight)



# Nearby genes for these loci?
all_sigmar_nearby_annotation %>%
  filter(marker %in% loci_highlight) %>%
  distinct(marker, nearby_annotation, within_annotations)



# # Allele frequencies for this marker in different populations
# allele_freq <- all_marker_freq %>%
#   filter(marker %in% loci_highlight)
#
# # Subset individuals that are homozygous for the minor/favorable allele
# favorable_genotype <- 1
# indivs_with_favorable_allele <- geno_mat_all[,loci_highlight] %>%
#   subset(., . == favorable_genotype)
#
# g <- pca_tidy %>%
#   left_join(., rownames_to_column(as.data.frame(geno_mat_all[,loci_highlight, drop = FALSE]), "individual")) %>%
#   rename_at(vars(all_of(loci_highlight)), ~paste0("marker", seq_along(.))) %>%
#   mutate_at(vars(contains("marker")), as.factor) %>%
#   ggplot(aes(x = PC1, y = PC2, color = marker1, shape = category, label = individual)) +
#   geom_point(size = 2) +
#   scale_x_continuous(name = paste0("PC1 (", round(pc_varexp["PC1"] * 100, 3), "%)"), breaks = pretty,
#                      labels = format_numbers) +
#   scale_y_continuous(name = paste0("PC2 (", round(pc_varexp["PC2"] * 100, 3), "%)"), breaks = pretty,
#                      labels = format_numbers) +
#   scale_alpha_discrete(guide = FALSE) +
#   theme_genetics(base_size = 12) +
#   theme(legend.position = "right", legend.justification = c(0.99, 0.99))
#
# plotly::ggplotly(g)












# Figure 5: soil texture association example --------------------------




## First show the locus on Chromosome 4 associated with multiple worldclim temperature
## variables

loci_highlight <- c("S02_38027635", "S03_4582157")

locus1_summary <- egwas_sigmar %>%
  filter(marker %in% loci_highlight) %>%
  # add environmental variable full names
  left_join(., select(eaa_environmental_vars, -class)) %>%
  arrange(marker)

# Only look at manhattans for topsoil and silt/sand/clay
locus1_summary1 <- locus1_summary %>%
  filter(str_detect(variable, "silt|sand|clay"), str_detect(variable, "topsoil"))

# Create a highlight window for this snp
locus1_highlight_rect <- snp_info %>%
  filter(marker %in% loci_highlight) %>%
  mutate(left = pos - 1000, right = pos + 1000)


# Create a merged manhattan plot
#
# Subset manhattan plots for biovars associated with this locus
manhattan_list <- subset(eaa_gwas_manhattan_plot_df, variable %in% locus1_summary1$variable, plot, drop = TRUE) %>%
  map(., ~. + geom_vline(data = locus1_highlight_rect, aes(xintercept = pos), lwd = 3, color = "grey75", alpha = 0.2) +
        theme(axis.title.y = element_blank())) %>%
  map(~{.x$labels$caption <- NULL; .x}) %>% # Remove the caption
  map(~{.x$layers[[2]]$aes_params$size <- 0.5; .x}) %>% # Modify point size
  # Modify theme
  map(~. + theme_genetics(base_size = base_font_size) +
        theme(panel.spacing.x = unit(0, "lines"), strip.placement = "outside",
              axis.text.x = element_blank(), axis.ticks.length.x = unit(0.25, "lines"), axis.line.x = element_blank(),
              legend.position = c(1, 1), legend.justification = c(1, 1), legend.direction = "vertical",
              legend.background = element_rect(fill = alpha("white", 0)),
              axis.title.y = element_blank(), strip.background.x = element_blank())) %>%
  # Modify the scale for linetype
  map(~ . + scale_linetype(name = NULL, guide = guide_legend(ncol = 1, keyheight = unit(0.5, "line")))) %>%
  modify_at(.x = ., .at = -length(.), ~. + theme(axis.ticks.x = element_blank(), strip.text.x = element_blank())) %>%
  modify_at(.x = ., .at = 1, ~ . + theme(legend.direction = "horizontal")) %>%
  modify_at(.x = ., .at = -1, ~. + theme(legend.position = "none"))

# Text grob
y_axis_label <- grid::textGrob(label = expression(-log[10](p-value)), rot = 90,
                               gp = grid::gpar(fontsize = base_font_size))

# Combine these
g_manhattan_combine <- cowplot::plot_grid(plotlist = manhattan_list, ncol = 1, align = "v",
                                          rel_heights = c(1, rep(0.7, length(manhattan_list)-2), 0.85))
g_manhattan_combine1 <- cowplot::plot_grid(y_axis_label, g_manhattan_combine, rel_widths = c(0.03, 1))

# Save
ggsave(filename = "soil_texture_ROI_manhattan_merge.jpg", plot = g_manhattan_combine1, path = fig_dir,
       height = 3.5, width = 3.5, dpi = 1000)


## Plot the geographic distribution of the alleles at these SNPs

allele_geo_list <- sigmar_allele_geo_dist %>%
  filter(marker %in% loci_highlight) %>%
  split(.$marker) %>%
  map("plot2_alt") %>%
  map(1)


## Bioclim variables versus marker alleles
bioclim_allele_list <- sigmar_summary %>%
  filter(marker %in% loci_highlight, str_detect(variable, "subsoil", negate = TRUE)) %>%
  split(.$marker) %>%
  map(~{
    plot_list <- .x$plot1
    # Combine plot data
    plot_data_combined <- map(plot_list, "data") %>%
      map2_dfr(.x = ., .y = map(plot_list, ~.x$labels$y), ~mutate(.x, variable = .y)) %>%
      mutate(genotypes_overall = genotype,
             variable = str_replace_all(variable, "\n", " "),
             variable = str_wrap(variable, width = 20))

    g_use <- plot_list[[1]]
    g_use$labels$y <- NULL
    g_use$labels$x <- NULL
    g_use$data <- plot_data_combined

    g_use1 <- g_use +
      facet_wrap(~ variable, nrow = 1, scales = "free_y") +
      theme(strip.background = element_blank())

    g_use1

  })


# Combine the geographic distribution plots and allele effect plots
sigmar_summary_soil_list <- map(names(allele_geo_list), ~{
  top_plot <- allele_geo_list[[.x]]
  bottom_plot <- bioclim_allele_list[[.x]]

  ref_allele <- subset(marker_metadata, marker == .x, ref_allele, drop = TRUE)
  alt_allele <- subset(marker_metadata, marker == .x, alt_allele, drop = TRUE)

  # Create genotype colors
  genotype_colors <- colorRampPalette(colors = allele_colors)(3) %>%
    setNames(., c(paste0(ref_allele, ref_allele), paste0(ref_allele, alt_allele), paste0(alt_allele, alt_allele)))

  # Edit subtitle
  top_plot1 <- top_plot +
    scale_color_manual(values = genotype_colors, name = paste0(.x, "\nmarker genotype"),
                       guide = guide_legend(override.aes = list(size = 2))) +
    theme_nothing(font_size = base_font_size) +
    theme(legend.position = c(0.75, 0.10), legend.justification = c(0,0))


  # Combine
  plot_grid(top_plot1, bottom_plot, ncol = 1, rel_heights = c(1, 0.8))

})

names(sigmar_summary_soil_list) <- names(allele_geo_list)

# Save plots
for (i in seq_along(sigmar_summary_soil_list)) {
  marker <- names(sigmar_summary_soil_list)[i]
  filename <- paste0("example_associations_marker_", marker, ".jpg")
  ggsave(filename = filename, plot = sigmar_summary_soil_list[[i]], path = fig_dir,
         height = 3, width = 3.5, dpi = 1000)

}


# Just save the geographic distribution and association with variables

# Merge using patchwork
combined_plot <- g_manhattan_combine1 + sigmar_summary_soil_list$S02_38027635 +
  plot_layout(widths = c(0.65, 1)) +
  plot_annotation(tag_levels = "A") +
  theme(plot.tag = element_text(size = base_font_size + 2))

# Save
ggsave2(filename = "figure5_soil_associations", plot = sigmar_summary_soil_list$S02_38027635,
       height = 3, width = 3.5, dpi = 500)


bioclim_allele_list$S02_38027635$data %>%
  group_by(variable, genotype) %>%
  summarize(biovar_mean = mean(biovar), nInd = n(), .groups = "drop") %>%
  filter(map_lgl(str_split(genotype, ""), ~n_distinct(.) == 1)) %>%
  mutate(genotype_class = ifelse(nInd == min(nInd), "minor", "major")) %>%
  select(-genotype, -nInd) %>%
  spread(genotype_class, biovar_mean) %>%
  mutate(diff = minor - major)

bioclim_allele_list$S03_4582157$data %>%
  group_by(variable, genotype) %>%
  summarize(biovar_mean = mean(biovar), nInd = n(), .groups = "drop") %>%
  filter(map_lgl(str_split(genotype, ""), ~n_distinct(.) == 1)) %>%
  mutate(genotype_class = ifelse(nInd == min(nInd), "minor", "major")) %>%
  select(-genotype, -nInd) %>%
  spread(genotype_class, biovar_mean) %>%
  mutate(diff = minor - major)



# Frequency of alleles in all populations
subset(all_marker_freq, marker %in% loci_highlight)


# Nearby genes for these loci?
all_sigmar_nearby_annotation %>%
  filter(marker %in% loci_highlight) %>%
  mutate(annotation_use = map2(within_annotations, nearby_annotation, ~{
    if(nrow(.x)>0) mutate(filter(.y, start == .x$start[1], end == .x$end[1]), min_distance = 0) else .y
  })) %>%
  distinct(marker, annotation_use) %>%
  unnest(annotation_use) %>%
  as.data.frame()


         annotation_use = modify_at(annotation_use)


           map2(within_annotations, nearby_annotation, ~{
    if (nrow(.x) > 0) mutate(.x, min_distance = 0) else as.data.frame(top_n(x = .y, n = 1, wt = -min_distance))
  })) %>%
  distinct(marker, annotation_use)







# What about the other variables
sigmar_summary %>%
  filter(marker %in% loci_highlight) %>%
  mutate(data = map(plot1, "data")) %>%
  select(variable, marker, full_name, data) %>%
  unnest(data) %>%
  group_by(variable, genotype) %>%
  summarize(biovar_mean = mean(biovar), nInd = n(), .groups = "drop") %>%
  filter(map_lgl(str_split(genotype, ""), ~n_distinct(.) == 1)) %>%
  mutate(genotype_class = ifelse(nInd == min(nInd), "minor", "major")) %>%
  select(-genotype, -nInd) %>%
  spread(genotype_class, biovar_mean) %>%
  mutate(diff = minor - major)


# Allele frequencies
subset(all_marker_freq, marker %in% loci_highlight)



# Nearby genes for these loci?
egwas_sigmar_nearby_annotation1 %>%
  filter(marker %in% loci_highlight) %>%
  select(-variable, -nearby_annotation) %>%
  distinct() %>%
  unnest(attributes) %>%
  as.data.frame()




# Supplemental Figure SXX - all chromosome LD -----------------------------

# Combine the full-length and 200 kbp cutoffs for all-chromosome LD

g_ld_decay_all_chromosomes <- g_ld_decay_all_chrom1 + g_ld_decay_all_chrom2 +
  plot_annotation(tag_levels = "A")


# Save
ggsave2(filename = "figureSXX_ld_decay_all_chromosomes", plot = g_ld_decay_all_chromosomes,
       width = 10, height = 4, dpi = 500)






# Supplemental figure SXX - correlation of environmental vars -------------

bioclim_dat_cor <- germplasm_bioclim_data %>%
  select(location_abbr, all_of(env_variables_of_interest)) %>%
  distinct() %>%
  as.data.frame() %>%
  column_to_rownames("location_abbr") %>%
  as.matrix() %>%
  cor()

upper_cor <- lower_cor <- bioclim_dat_cor
upper_cor[!upper.tri(upper_cor, diag = T)] <- NA
lower_cor[!lower.tri(lower_cor, diag = T)] <- NA

heatmap_dat <- heatmap(bioclim_dat_cor)
variable_levels <- row.names(bioclim_dat_cor)[heatmap_dat$rowInd]

g_cor <- GGally::ggcorr(data = bioclim_dat_cor[variable_levels, variable_levels])

g_cor1 <- bind_rows(mutate(g_cor$data, upper = FALSE), mutate(rename(g_cor$data, y = x, x = y), upper = TRUE)) %>%
  add_row(x = variable_levels, y = variable_levels, coefficient = 1, label = NA) %>%
  mutate_at(vars(x, y), ~factor(., levels = variable_levels)) %>%
  mutate(coefficient = ifelse(upper, NA, coefficient),
         label = ifelse(upper, label, NA)) %>%
  ggplot(aes(x = x, y = y, fill = coefficient, label = label)) +
  geom_tile() +
  geom_text(size = 1.5) +
  scale_fill_gradient2(low = scales::muted("blue"), high = scales::muted("red"), na.value = "white",
                       limits = c(-1, 1), name = "Correlation") +
  theme_genetics() +
  theme(axis.title = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))

# Save
ggsave2(filename = "figureSXX_biovar_correlation", plot = g_cor1, height = 5, width = 7.5, dpi = 500)



# Supplemental Figure SXX: dissection of Fst outliers ---------------------

# Annotation for Fst outliers
fst_outliers_nearby_annotation1 %>%
  select(-nearby_annotation) %>%
  unnest(attributes) %>%
  as.data.frame()

loci_highlight <- fst_outliers$marker


allele_geo_list <- sigmar_allele_geo_dist %>%
  filter(marker %in% loci_highlight) %>%
  split(.$marker) %>%
  map("plot2_alt") %>%
  map(1)

allele_freq <- sigmar_allele_geo_dist %>%
  filter(marker %in% loci_highlight) %>%
  select(marker, allele_freq)


## Bioclim variables versus marker alleles
bioclim_allele_list <- sigmar_summary %>%
  filter(marker %in% loci_highlight, !is.na(full_name)) %>%
  split(.$marker) %>%
  map(~{
    plot_list <- .x$plot1
    # Combine plot data
    plot_data_combined <- map(plot_list, "data") %>%
      map2_dfr(.x = ., .y = map(plot_list, ~.x$labels$y), ~mutate(.x, variable = .y)) %>%
      mutate(genotypes_overall = genotype,
             variable = str_replace_all(variable, "\n", " "),
             variable = str_wrap(variable, width = 20))

    g_use <- plot_list[[1]]
    g_use$labels$y <- NULL
    g_use$labels$x <- NULL
    g_use$data <- plot_data_combined

    g_use1 <- g_use +
      facet_wrap(~ variable, nrow = 1, scales = "free_y") +
      theme(strip.background = element_blank())

    g_use1

  })




# Combine the geographic distribution plots and allele effect plots
sigmar_summary_list <- map(names(allele_geo_list), ~{
  top_plot <- allele_geo_list[[.x]]
  bottom_plot <- bioclim_allele_list[[.x]]

  allele_freq_i <- unnest(subset(allele_freq, marker == .x))
  ref_allele <- as.character(subset(allele_freq_i, allele_type == "ref", allele, drop = TRUE)[1])
  alt_allele <- as.character(subset(allele_freq_i, allele_type == "alt", allele, drop = TRUE)[1])

  # Create genotype colors
  genotype_colors <- colorRampPalette(colors = allele_colors)(3) %>%
    setNames(., c(paste0(ref_allele, ref_allele), paste0(ref_allele, alt_allele), paste0(alt_allele, alt_allele)))

  # Edit subtitle
  top_plot1 <- top_plot +
    scale_color_manual(values = genotype_colors, name = paste0(.x, "\nmarker genotype"),
                       guide = guide_legend(override.aes = list(size = 2))) +
    theme_nothing(font_size = base_font_size) +
    theme(legend.position = c(0.75, 0.10), legend.justification = c(0,0))


  # Combine
  plot_grid(top_plot1, bottom_plot, ncol = 1, rel_heights = c(1, 0.8)) +
    theme(panel.background = element_rect())

})

names(sigmar_summary_list) <- names(allele_geo_list)

# Combine all figures
g_combined <- plot_grid(plotlist = sigmar_summary_list, ncol = 2, labels = LETTERS[seq_along(sigmar_summary_list)])

ggsave2(filename = "figureSXX_Fst_outliers", plot = g_combined, height = 10, width = 8, dpi = 500)




# Supplemental Figure SXX: dissection of SPA outliers ---------------------


loci_highlight <- spa_outliers$marker


allele_geo_list <- sigmar_allele_geo_dist %>%
  filter(marker %in% loci_highlight) %>%
  split(.$marker) %>%
  map("plot2_alt") %>%
  map(1)

allele_freq <- sigmar_allele_geo_dist %>%
  filter(marker %in% loci_highlight) %>%
  select(marker, allele_freq)


## Bioclim variables versus marker alleles
bioclim_allele_list <- sigmar_summary %>%
  filter(marker %in% loci_highlight, !is.na(full_name)) %>%
  split(.$marker) %>%
  map(~{
    plot_list <- .x$plot1
    # Combine plot data
    plot_data_combined <- map(plot_list, "data") %>%
      map2_dfr(.x = ., .y = map(plot_list, ~.x$labels$y), ~mutate(.x, variable = .y)) %>%
      mutate(genotypes_overall = genotype,
             variable = str_replace_all(variable, "\n", " "),
             variable = str_wrap(variable, width = 20))

    g_use <- plot_list[[1]]
    g_use$labels$y <- NULL
    g_use$labels$x <- NULL
    g_use$data <- plot_data_combined

    g_use1 <- g_use +
      facet_wrap(~ variable, nrow = 1, scales = "free_y") +
      theme(strip.background = element_blank())

    g_use1

  })




# Combine the geographic distribution plots and allele effect plots
sigmar_summary_list <- map(names(allele_geo_list), ~{
  top_plot <- allele_geo_list[[.x]]
  bottom_plot <- bioclim_allele_list[[.x]]

  allele_freq_i <- unnest(subset(allele_freq, marker == .x))
  ref_allele <- as.character(subset(allele_freq_i, allele_type == "ref", allele, drop = TRUE)[1])
  alt_allele <- as.character(subset(allele_freq_i, allele_type == "alt", allele, drop = TRUE)[1])

  # Create genotype colors
  genotype_colors <- colorRampPalette(colors = allele_colors)(3) %>%
    setNames(., c(paste0(ref_allele, ref_allele), paste0(ref_allele, alt_allele), paste0(alt_allele, alt_allele)))

  # Edit subtitle
  top_plot1 <- top_plot +
    scale_color_manual(values = genotype_colors, name = paste0(.x, "\nmarker genotype"),
                       guide = guide_legend(override.aes = list(size = 2))) +
    theme_nothing(font_size = base_font_size) +
    theme(legend.position = c(0.75, 0.10), legend.justification = c(0,0))


  # Combine
  plot_grid(top_plot1, bottom_plot, ncol = 1, rel_heights = c(1, 0.8)) +
    theme(panel.background = element_rect())

})

names(sigmar_summary_list) <- names(allele_geo_list)

# Combine all figures
g_combined <- plot_grid(plotlist = sigmar_summary_list, ncol = 2, labels = LETTERS[seq_along(sigmar_summary_list)])

ggsave2(filename = "figureSXX_SPA_outliers", plot = g_combined, height = 10, width = 8, dpi = 500)










# Supplemental Table S01: Name and description of environmental variables ---------------

eaa_environmental_vars %>%
  arrange(class, variable) %>%
  rename_all(~str_replace_all(., "_", " ") %>% str_to_title()) %>%
  write_csv(x = ., file = file.path(fig_dir, "table_S01_env_variable_description.csv"))




# Supplemental Table S02: all EAA significant markers ---------------------------------


# Add the full names of the environmental variables
egwas_sigmar %>%
  select(-model, -contains("score")) %>%
  left_join(., eaa_environmental_vars) %>%
  select(variable, class, marker, chrom, pos, p_value) %>%
  rename_all(~str_replace_all(., "_", " ") %>% str_to_title()) %>%
  # Save as a CSV
  write_csv(x = ., file = file.path(fig_dir, "table_S02_egwas_significant_associations.csv"))



# Supplemental Table S03: SPA and Fst outliers ---------------------------------

all_sigmar %>%
  filter(! variable %in% eaa_environmental_vars$variable) %>%
  select(variable, marker, chrom, pos, score) %>%
  rename_all(~str_replace_all(., "_", " ") %>% str_to_title()) %>%
  # Save as a CSV
  write_csv(x = ., file = file.path(fig_dir, "table_S03_spa_fst_outlier_markers.csv"))









# Code for other results described in the paper ---------------------------

germ_meta1 <- germ_meta %>%
  left_join(., distinct(eaa_environmental_data, origin_name, location_abbr))

## Number of individuals per location
aggregate(individual ~ location_abbr, data = germ_meta1, FUN = n_distinct) %>%
  arrange(individual)



# Summarize within-location pairwise distances
pairwise_dist_data %>%
  filter(location1 == location2) %>%
  group_by(location1) %>%
  summarize(mean_dist = mean(distance), min_dist = min(distance), max_dist = max(distance),
            range_dist = max_dist - min_dist) %>%
  left_join(., aggregate(individual ~ location_abbr, data = germ_meta1, FUN = n_distinct), by = c("location1" = "location_abbr")) %>%
  arrange(desc(mean_dist)) %>%
  mutate(mean_similarity = 1 - mean_dist)


## eGWAS results ##

# Subset the monthly variables
egwas_sigmar1 <- egwas_sigmar %>%
  filter(str_detect(variable, "prec|tmax|tmin|tmean", negate = TRUE))

# Number of total marker associations
nrow(egwas_sigmar1)

# Number of unique variables and markers
n_distinct(egwas_sigmar1$variable)
n_distinct(egwas_sigmar1$marker)


# Determine significant

# Number of markers per variable
xtabs(~ variable, egwas_sigmar1) %>%
  sort()

# Number of variables per marker
xtabs(~ marker, egwas_sigmar1) %>%
  sort()

# Breakdown by variable type
egwas_sigmar1 %>%
  distinct(variable) %>%
  left_join(., eaa_environmental_vars) %>%
  xtabs(~ class, .)

xtabs(~ class, eaa_environmental_vars)


## Overlap of Fst and egwas snps
eaa_marker_fst <- global_marker_fst %>%
  mutate(Fst_rank = rank(-Fst)) %>%
  filter(marker %in% unique(egwas_sigmar$marker)) %>%
  arrange(desc(Fst)) %>%
  as.data.frame()

mean(global_marker_fst$Fst)
sd(global_marker_fst$Fst)

mean(eaa_marker_fst$Fst)
sd(eaa_marker_fst$Fst)

# Compare distributions
(mw_test_out <- wilcox.test(x = eaa_marker_fst$Fst, y = global_marker_fst$Fst, conf.int = T))

## Highlight a case example

## First show the locus on Chromosome 4 associated with multiple worldclim temperature
## variables

loci_highlight <- eaa_marker_fst$marker[1:4]

egwas_sigmar %>%
  filter(marker %in% loci_highlight) %>%
  # add environmental variable full names
  left_join(., select(eaa_environmental_vars, -class))


## Bioclim variables versus marker alleles
bioclim_allele_list <- sigmar_summary %>%
  filter(marker %in% loci_highlight) %>%
  split(.$marker) %>%
  map(~{
    plot_list <- .x$plot1
    # Combine plot data
    plot_data_combined <- map(plot_list, "data") %>%
      map2_dfr(.x = ., .y = map(plot_list, ~.x$labels$y), ~mutate(.x, variable = .y)) %>%
      mutate(genotypes_overall = genotype,
             variable = str_replace_all(variable, "\n", " "),
             variable = str_wrap(variable, width = 20))

    g_use <- plot_list[[1]]
    g_use$labels$y <- NULL
    g_use$labels$x <- NULL
    g_use$data <- plot_data_combined

    g_use1 <- g_use +
      facet_wrap(~ variable, nrow = 1, scales = "free_y") +
      theme(strip.background = element_blank())

    g_use1

  })







# Highlight results for two environmental variables -----------------------

# # Load the bioclim rasters
# raster_list <- list.files(path = file.path(data_dir, "rasterFiles"), full.names = TRUE) %>%
#   lapply(raster)
#
# names(raster_list) <- sapply(raster_list, names) %>%
#   str_remove(., "_data_raster") %>%
#   str_remove(., "soilgrids_|worldclim_")
#
# # Aggregate elevation and soilgrids rasters
# raster_list_aggregates <- raster_list %>%
#   modify_at(., subset(eaa_environmental_vars, class == "soil", variable, drop = TRUE),
#             ~aggregate(x = ., fact = floor(res(raster_list$bio1) / res(raster_list$bdod_subsoil)), mean))
#
#
# bioclim_raster_df_list <- lapply(X = raster_list, FUN = as.data.frame, xy = TRUE) %>%
#   map(as_tibble) %>%
#   imap(~`names<-`(.x, c("x", "y", .y)))
#
#
#
#
# # Highlight example results for presentations -----------------------------
#
# biovar_select <- c("bio1", "phh2o_topsoil", "sand_topsoil")
#
#
#
#
# ##
# ## Second look at soil pH and sand content associations
# ##
# loci_summary <- sigmar_summary %>%
#   filter(variable %in% biovar_select,
#          str_detect(variable, "bio", negate = TRUE),
#          class == "mlmm") %>%
#   left_join(., mlmm_out_list[biovar_select] %>% map(1) %>% imap_dfr(~mutate(.x, variable = .y))) %>%
#   group_by(variable) %>%
#   top_n(x = ., n = 1, wt = -p_value) %>%
#   # # Pull the fdr thresholds
#   # left_join(., mutate(eaa_gwas_manhattan_plot_df, fdr = map_dbl(plot, ~min(.x$layers[[1]]$data$fdr_score, na.rm = TRUE)))) %>%
#   # # Highlight regions above -log10(p) >= 4
#   # filter(-log10(p_value) >= fdr) %>%
#   # Manual filter
#   filter(
#     (variable == "sand_topsoil" & str_detect(term, "S02|S03|S06")) | variable == "phh2o_topsoil",
#     p_value < 0.0025
#   )
#
#
# # Subset the mahattan plots
# #
# #
# # Create highlights of the loci
# loci2_highlight_rect <- snp_info %>%
#   inner_join(., select(loci_summary, variable, term), by = c("marker" = "term")) %>%
#   mutate(left = pos - 10000, right = pos + 10000) %>%
#   group_by(variable) %>%
#   nest() %>%
#   ungroup()
#
# manhattan_list <- eaa_gwas_manhattan_plot_df %>%
#   inner_join(., loci2_highlight_rect) %>%
#   mutate(plot = map2(.x = plot, .y = data, ~{
#     .x1 <- .x + geom_vline(data = .y, aes(xintercept = pos), lwd = 5, color = "grey75", alpha = 0.2) +
#       theme(axis.title.y = element_blank(), legend.direction = "horizontal")
#     .x1$labels$subtitle <- str_replace(string = .x1$labels$subtitle, pattern = "Ph", replacement = "pH")
#     .x1$layers[[2]]$aes_params$size <- 1.5 # Modify point size
#     .x1
#   }))
#
# # Save these
# for (i in seq_len(nrow(manhattan_list))) {
#   filename <- paste0(manhattan_list$variable[i], "_manhattan_toplot.jpg")
#   ggsave(filename = filename, plot = manhattan_list$plot[[i]], path = fig_dir,
#          width = 10, height = 3, dpi = 1000)
# }
#
#
# # Subset the effect plots
# effect_plot_list <- setNames(object = loci_summary$plot1, nm = loci_summary$variable) %>%
#   map(., ~. + theme_genetics(base_size = 16) + theme(axis.title.x = element_blank())) %>%
#   map(., ~. + ylab(str_wrap(str_trim(str_split(string = .x$labels$y, pattern = "\\(")[[1]][1]), 20))) %>%
#   map(., ~{
#     .x$labels$y <- str_replace(string = .x$labels$y, pattern = "Ph", replacement = "pH")
#     .x$layers[[2]]$aes_params$size <- 1
#     .x
#   })
#
# # Subset the map plots
# map_plot_list <- setNames(object = loci_summary$plot2, nm = loci_summary$variable) %>%
#   map(., ~{
#     .x + theme_void(base_size = 16) +
#       # labs(subtitle = "Spatial allele frequencies") +
#       theme(legend.position = c(0.90, 0.05), legend.justification = c(1,0))
#   })
#
# # Merge plots
# loci_effect_map_merge <- map2(.x = map_plot_list, .y = effect_plot_list, ~{
#   cowplot::plot_grid(.x, .y, rel_widths = c(1, 0.6))
# })
#
# # Save these
# for (i in seq_along(loci_effect_map_merge)) {
#   filename <- paste0(names(loci_effect_map_merge)[i], "_map_effect_toplot.jpg")
#   ggsave(filename = filename, plot = loci_effect_map_merge[[i]], path = fig_dir,
#          width = 10, height = 2, dpi = 1000)
# }
#
#
#
#
# # Other summary statistics
# all_assoc <- filter(sigmar_summary, class == "mlmm")
# nrow(all_assoc)
# n_distinct(all_assoc$variable)
# n_distinct(all_assoc$term)
#
# xtabs(~ chrom, snp_info, subset = marker %in% all_assoc$term)
#
