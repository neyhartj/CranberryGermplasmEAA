# CranberryGermplasmEAA
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
library(khroma)


# Load the startup script
source("startup.R")

# Labels for subfigures
subfigure_labels <- LETTERS

# Base font size
base_font_size <- 8

# Settings for figures
fig_dpi <- 300
fig_device <- "tiff"


# Read in data ------------------------------------------------------------

# Read in the base data
load(file.path(data_dir, "population_metadata_and_genotypes.RData"))
load(file.path(data_dir, "germplasm_origin_bioclim_data.RData"))
load(file.path(result_dir, "population_genetics_stats.RData"))
# Load the eGWAS results
load(file.path(result_dir, "eaa_gwas_results.RData"))
load(file.path(result_dir, "egwas_sigmar_analysis.RData"))


# Load the annotation data
gff_file <- file.path(data_dir, "Vaccinium_macrocarpon_BenLear_v2_annotations.gff")
# read in the GFF file
cranberry_gff <- ape::read.gff(file = gff_file)
cran_gene_ranges <- rtracklayer::import.gff3(con = gff_file)


# Calculate min/max chrom lengths
chromlengths <- cranberry_gff %>%
  filter(type == "chromosome", seqid != "Unknown") %>%
  select(chrom = seqid, start, end) %>%
  mutate(chrom = str_pad(chrom, 2, pad = "0")) %>%
  arrange(chrom)

# Rename genotype matrices
geno_mat_wild <- geno_mat_wild_filter3
geno_hmp_wild <- geno_hmp_wild_filter3
geno_mat_all <- geno_mat_all_filter3



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
  summarize(nInd = n(), ref_allele_freq = mean(ref_allele_count, na.rm = TRUE) / 2, .groups = "drop")


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
#   ggsave(filename = filename, plot = combined_plot1, path = fig_dir, width = 8, height = 2.5, dpi = 150)
#
#   pb$tick()
#
# }





# Figure 1: geographic origin -------------

# Plot origin of wild germplasm

wild_germ_meta1 <- pop_metadata %>%
  filter(category == "Wild")


# Create a box for x and y axes
xlim <- range(pretty(wild_germ_meta1$longitude))
ylim <- range(pretty(wild_germ_meta1$latitude))

# Breaks for population size
pop_size_breaks <- pretty(range(wild_germ_meta_origin_summ$nSamples), n = 3)

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
  # labs(subtitle = "Wild cranberry collection sites") +
  theme_void(base_size = 14) +
  theme(legend.position = c(0.90, 0.05), legend.justification = c(1,0), legend.box = "vertical", legend.box.just = "left",
        panel.border = element_rect(fill = NA))

# Save the figure
ggsave2(filename = "figure1_cranberry_wild_germplasm_origins1", plot = g_map, path = fig_dir,
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

# Plot variance explained
plot(K_pca)

# Cumulative sum
plot(cumsum(pc_varexp))
head(pc_varexp)
sum(pc_varexp[1:3])

# Calculate the median position of each color and plot a label for that
group_label <- pca_tidy %>%
  group_by(state) %>%
  summarize_at(vars(contains("PC")), median)

# Pairwise plots of PCs 1-3
pc_pairs <- combn(x = paste0("PC", 1:3), m = 2, simplify = FALSE)
names(pc_pairs) <- sapply(pc_pairs, paste0, collapse = "_")
pairwise_pc_plots <- map(pc_pairs, ~{
  x <- .x[1]
  y <- .x[2]

  ggplot(data = pca_tidy, aes_string(x = x, y = y, color = "state")) +
    geom_point(aes(alpha = location_abbr, linetype = individual)) +
    ggrepel::geom_label_repel(data = group_label, aes(label = state), fill = NA, key_glyph = "point", size = 2,
                              max.overlaps = 15, point.padding = unit(10, "line"), box.padding = unit(0.2, "line")) +
    scale_color_discrete(guide = FALSE) +
    scale_x_continuous(name = paste0(x, " (", round(pc_varexp[x] * 100, 3), "%)"), breaks = pretty,
                       labels = format_numbers) +
    scale_y_continuous(name = paste0(y, " (", round(pc_varexp[y] * 100, 3), "%)"), breaks = pretty,
                       labels = format_numbers) +
    scale_alpha_discrete(guide = FALSE) +
    theme_genetics() +
    theme(legend.position = c(0.99, 0.99), legend.justification = c(0.99, 0.99))

})


# add Origin and plot
pc_plot <- pairwise_pc_plots$PC1_PC2

# Save this
ggsave(filename = "snp_pca_wild_germplasm.jpg", plot = pc_plot, path = fig_dir, width = 5, height = 4, dpi = 1000)




## Plot PCA of all entries

# Gather the eigenvector
pca_tidy_all <- broom::tidy(K_pca_all) %>%
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
group_label <- pca_tidy_all %>%
  group_by(state) %>%
  summarize_at(vars(contains("PC")), median)

# add Origin and plot
pc_plot_all <- pca_tidy_all %>%
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

# Plot decay curves of all chromosomes

# Unfiltered distance
g_ld_decay_curves_all_chrom <- ld_wild_decay1 %>%
  ggplot(aes(x = d / 1000, y = yhat, color = LG)) +
  geom_line(lwd = 0.75) +
  scale_x_continuous(name = "Distance (kb)", breaks = pretty) +
  scale_y_continuous(name = expression("Linkage disequilibrium ("*r^2*")")) +
  scale_color_discreterainbow(name = "Chrom.", guide = guide_legend(ncol = 2, keyheight = unit(0.5, "line"))) +
  # facet_wrap(~ LG, labeller = labeller(LG = function(x) paste0("Chrom. ", x))) +
  theme_genetics(base_size = base_font_size) +
  theme(legend.position = c(0.35, 0.8))

# Filter distance
g_ld_decay_curves_all_chrom_filtered <- g_ld_decay_curves_all_chrom %>%
  modify_at("data", ~filter(., d <= 2e5))

# Combine plots

# Determine the left and right plots
left_plot <- g_ld_decay_curves_all_chrom_filtered
right_plot <- g_ld_decay_curves_all_chrom +
  theme(legend.position = "none", axis.title = element_blank(),
        plot.background = element_rect(fill = alpha("white", 1), color = "black"))

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

# Minimum BP of separation before segments are labeled
min_bp_dist <- 1e4

# Assign class labels to SNPs within a certain distance
all_stats_to_plot_unique_labels <- all_stats_to_plot_unique %>%
  arrange(chrom, pos) %>%
  group_by(chrom) %>%
  do({
    df <- .
    df1 <- mutate(df, pos_diff = diff(c(0, pos)))
    which_label <- union(which(df1$pos_diff < min_bp_dist), which(df1$pos_diff < min_bp_dist)-1)
    df1$label = as.character(NA)
    df1$label[which_label] <- as.character(df1$class[which_label])
    df1
  }) %>% ungroup() %>%
  mutate(nudge_y = rep(c(4, -4), length.out = nrow(.)), nudge_y = ifelse(is.na(label), as.numeric(NA), nudge_y))

# Segment width adjuster
segWidth <- 200000

bp_div <- 1e6


# A transposed version
g_chromlengths_rotate <- chromlengths1 %>%
  ggplot(aes(x = start / bp_div, xend = end / bp_div, y = as.factor(0), yend = as.factor(0))) +
  geom_segment(lwd = 5, lineend = "round", color = "grey85") +
  # Add annotations
  geom_segment(data = subset(all_stats_to_plot_unique),
               aes(x = (pos - segWidth) / bp_div, xend = (pos + segWidth) / bp_div, color  = class,
                   alpha = as.factor(nVar)), lwd = 5) +
  geom_text_repel(data = all_stats_to_plot_unique_labels, aes(y = as.factor(0), x = pos / bp_div, label = label, color = class),
                  inherit.aes = FALSE, direction = "x", size = 2, hjust = 0, parse = TRUE, show.legend = FALSE,
                  nudge_y = all_stats_to_plot_unique_labels$nudge_y, force_pull = 0, force = 0.25,
                  min.segment.length = 0, ) +
  facet_grid(chrom ~ ., switch = "y", labeller = labeller(chrom = as_mapper(~paste("Chrom.", as.numeric(.))))) +
  scale_alpha_manual(guide = FALSE, values = seq(0.5, 1, length.out = 6)) +
  scale_color_muted(name = NULL, guide = guide_legend(nrow = 1, keywidth = unit(0, "line"), override.aes = list(lwd = 2)),
                    labels = function(x) parse(text = x)) +
  scale_y_discrete(name = NULL, labels = NULL) +
  scale_x_continuous(name = "Position (Mb)", breaks = pretty, expand = expansion(0.05, 0.15)) +
  theme_genetics(base_size = base_font_size) +
  theme(strip.text.y.left = element_text(size = base_font_size, angle = 0), strip.background.y = element_blank(),
        panel.spacing = unit(0, "line"), legend.position = "bottom", axis.line.y = element_blank(),
        axis.ticks.y = element_blank(), legend.box.margin = margin(), legend.margin = margin())

ggsave2(filename = "figure3_genomic_statistics_chroms_rotated", plot = g_chromlengths_rotate,
        width = 4, height = 5.5, dpi = 500)





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
  scale_linetype_discrete(name = NULL, guide = guide_legend(ncol = 1, keyheight = unit(0.5, "lines"))) +
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
  scale_x_continuous(name = "Chromosome 4 position (kb)", breaks = pretty) +
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
      ggplot(aes(x = genotypes_overall, y = biovar, group = genotypes_overall)) +
      geom_jitter(aes(fill = genotypes_overall), pch = 21, size = 1,
                  position = position_jitterdodge(jitter.width = 0.10, dodge.width = 0.75)) +
      geom_point(data = genotype_means, position = position_dodge(0.75), pch = "+", size = 5, color = neyhart_palette()[3]) +
      facet_wrap(~ variable, ncol = 3, scales = "free",
                 labeller = labeller(variable = as_mapper(~f_rename_variables(., units = T)))) + # dir = "v" means fill by column
      scale_fill_manual(values = genotype_colors, guide = FALSE) +
      scale_y_continuous(name = NULL, breaks = pretty) +
      scale_x_discrete(name = NULL) +
      theme_grey(base_size = base_font_size) +
      theme(strip.background = element_blank(), panel.grid.minor.x = element_blank(),
            panel.grid.major.x = element_blank(), legend.position = c(0.80, 0.25))

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


# Frequency of alleles in all populations
subset(all_marker_freq, marker %in% loci_highlight)



# Nearby genes for these loci?
all_sigmar_nearby_annotation %>%
  filter(marker %in% loci_highlight) %>%
  distinct(marker, nearby_annotation, within_annotations)



# Figure 5: soil texture association example --------------------------




## First show the locus on Chromosome 4 associated with multiple worldclim temperature
## variables

loci_highlight <- c("S02_38027635", "S03_4582157")
locus_highlight <- "S03_4582157"



egwas_sigmar %>%
  filter(marker %in% locus_highlight) %>%
  # add environmental variable full names
  left_join(., select(eaa_environmental_vars, -class))

locus1_summary <- egwas_sigmar %>%
  filter(marker %in% locus_highlight, str_detect(variable, "subsoil", negate = TRUE)) %>%
  # add environmental variable full names
  left_join(., select(eaa_environmental_vars, -class))

# Create a highlight window for this snp
locus1_highlight_rect <- snp_info %>%
  filter(marker %in% locus_highlight) %>%
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
  scale_linetype_discrete(name = NULL, guide = guide_legend(ncol = 1, keyheight = unit(0.5, "lines"))) +
  theme_genetics(base_size = base_font_size) +
  theme(panel.spacing.x = unit(0, "lines"), strip.placement = "outside",
        axis.text.x = element_blank(), axis.ticks.length.x = unit(0.25, "lines"), axis.line.x = element_blank(),
        legend.position = c(0.99, 0.99), legend.justification = c(0.99, 0.99), legend.direction = "horizontal",
        legend.background = element_rect(fill = alpha("white", 0)),
        strip.background.x = element_blank())


## Plot Part B - Genomic range of the gene

loci_grange <- subset(snp_info, marker %in% locus_highlight) %>%
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
  summarize(x = (max(x) + min(x))/2, y = min(y) - (2*offset), strand = first(strand)) %>%
  # mutate(x = modify_at(x, 3, ~. + 20000)) %>%
  mutate(y1 = as.numeric(as.factor(strand)))

# Nudge x factor
n_x <- rep(0, nrow(loci_overlaps_labels))
n_x[3] <- 10


g_loci_overlaps <- loci_overlaps_polygons %>%
  ggplot(aes(x = x / 1e3, y = y)) +
  geom_polygon(aes(fill = strand, group = group)) +
  geom_vline(data = locus1_highlight_rect, aes(xintercept = pos / 1e3), lwd = 2.5, color = "grey75", alpha = 0.2) +
  # geom_text(data = loci_overlaps_labels, mapping = aes(label = Name), size = 2, hjust = "inward") +
  geom_text_repel(data = loci_overlaps_labels, mapping = aes(label = Name, y = y1), size = 2, direction = "x",
                  hjust = "inward", nudge_y = ifelse(loci_overlaps_labels$strand == "+", 0.2, -0.2), nudge_x = n_x,
                  box.padding = unit(0.8, "line"), force_pull = 0, force = 1, min.segment.length = unit(0.5, "line")) +
  scale_fill_discrete(guide = FALSE) +
  # Adjust 'expand' argument to align the vertical lines from the manhattan plot and the gene model plot
  scale_x_continuous(name = "Chromosome 3 position (kb)", breaks = pretty) +
  scale_y_continuous(expand = expansion(mult = 0.1)) +
  theme_genetics(base_size = base_font_size) +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.line.y = element_blank())



## Part C - geographic distribution of the marker alleles

# Subset this from the list of plots above
allele_geo_list <- sigmar_summary %>%
  filter(marker %in% locus_highlight) %>%
  split(.$marker) %>%
  map("plot2") %>%
  map(1) %>%
  map(~.x + theme_nothing(font_size = base_font_size) + theme(legend.position = .x$theme$legend.position))


## Part D - Bioclim variables versus marker alleles

bioclim_allele_list <- sigmar_summary %>%
  filter(marker %in% loci_highlight) %>%
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
      ggplot(aes(x = genotypes_overall, y = biovar, group = genotypes_overall)) +
      geom_jitter(aes(fill = genotypes_overall), pch = 21, size = 1,
                  position = position_jitterdodge(jitter.width = 0.10, dodge.width = 0.75)) +
      geom_point(data = genotype_means, position = position_dodge(0.75), pch = "+", size = 5, color = neyhart_palette()[3]) +
      facet_wrap(~ variable, ncol = 3, scales = "free",
                 labeller = labeller(variable = as_mapper(~f_rename_variables(., units = T)))) + # dir = "v" means fill by column
      scale_fill_manual(values = genotype_colors, guide = FALSE) +
      scale_y_continuous(name = NULL, breaks = pretty) +
      scale_x_discrete(name = NULL) +
      theme_grey(base_size = base_font_size) +
      theme(strip.background = element_blank(), panel.grid.minor.x = element_blank(),
            panel.grid.major.x = element_blank(), legend.position = c(0.80, 0.25))

  })


## Combine everything

p1 <- g_man_combined
p2 <- g_loci_overlaps
p3 <- allele_geo_list$S03_4582157
p4 <- bioclim_allele_list$S03_4582157

# For patchwork, you need to specify the nesting design, as demonstrated below
g_combined <- ( (p1 + p2 + plot_layout(ncol = 1, heights = c(1, 0.12))) |
                  (p3 + p4 + plot_layout(ncol = 1, heights = c(1, 1))) ) +
  plot_layout(ncol = 2, widths = c(0.40, 1)) +
  plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(size = base_font_size + 2, face = "plain"))


# Save
ggsave2(filename = "figure5_chrom3_soil_precip_associations", plot = g_combined,
        width = 8, height = 6, dpi = 500)




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



# Allele frequencies
subset(all_marker_freq, marker %in% loci_highlight)



# Table 1: pairwise genetic distance -----------------------

# Summarize within-location pairwise distances
table1_prep <- pairwise_dist_data %>%
  filter(location1 == location2) %>%
  group_by(location1) %>%
  summarize(mean_distance = mean(distance), min_dist = min(distance), max_distance = max(distance)) %>%
  left_join(., aggregate(individual ~ location_abbr, data = germplasm_bioclim_data, FUN = n_distinct), by = c("location1" = "location_abbr")) %>%
  arrange(desc(mean_distance)) %>%
  select(location = location1, sample_size = individual, mean_pairwise_distance = mean_distance,
         min_pairwise_distance = min_dist, max_pairwise_distance = max_distance)


# Format for the manuscript
# Add location full names and state information
table1_prep %>%
  mutate_at(vars(contains("distance")), format_numbers) %>%
  mutate(annotation = paste0(mean_pairwise_distance, " (", min_pairwise_distance, ", ", max_pairwise_distance, ")")) %>%
  left_join(., mutate(distinct(eaa_environmental_data, location_abbr, location, state), location = str_to_title(location)), by = c("location" = "location_abbr")) %>%
  select(location = location.y, state, abbreviation = location, N = sample_size, annotation) %>%
  rename_all(~str_replace_all(., "_", " ") %>% str_to_title()) %>%
  # Save as a CSV
  write_csv(x = ., file = file.path(fig_dir, "table1_pairwise_genetic_distance.csv"))


# Supplemental Figure SXX: other PCA plots -------------------------------

# Plot combinations of PCs as a supplemental figure
g_combined_pc_plots <- wrap_plots(pairwise_pc_plots, ncol = 1)

# Save
ggsave(filename = "figureSXX_population_structure_pairwise.jpg", plot = g_combined_pc_plots,
       path = fig_dir, height = 8, width = 3, dpi = 500)



# Supplemental Figure S01: boxplots of environmental variables -----------

bioclim_data_plot <- germplasm_bioclim_data %>%
  as_tibble() %>%
  select(-individual) %>%
  distinct() %>%
  gather(variable, value, all_of(env_variables_of_interest)) %>%
  left_join(., eaa_environmental_vars)

# List of plots
bioclim_plot_list <- bioclim_data_plot %>%
  mutate(variable = factor(variable, levels = env_variables_of_interest)) %>%
  arrange(variable) %>%
  split(.$variable) %>%
  map(~{
    ggplot(data = .x, aes(x = variable, y = value, color = location_abbr)) +
      geom_boxplot(aes(group = 1), width = 0.5, outlier.shape = NA) +
      geom_jitter(width = 0.1) +
      geom_text_repel(aes(label = location_abbr), nudge_x = rep(c(-0.5, 0.5), length.out = nrow(.x)), hjust = 1, direction = "y",
                      min.segment.length = unit(100, "line"), size = 2) +
      scale_color_discrete(guide = FALSE) +
      scale_y_continuous(name = parse(text = str_replace_all(unique(.x$unit), " ", "~")), breaks = pretty) +
      scale_x_discrete(name = NULL, labels = f_rename_variables) +
      theme_genetics(base_size = base_font_size)
  })

# Combine
g_bioclim_boxplots <- wrap_plots(bioclim_plot_list, ncol = 5)

# Save
ggsave(filename = "figureS01_bioclim_boxplots.jpg", plot = g_bioclim_boxplots,
       path = fig_dir, height = 25, width = 10, dpi = 500)


# Supplemental Figure S02: correlation of environmental vars -------------

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
ggsave(filename = "figureS02_biovar_correlation.jpg", plot = g_cor1, path = fig_dir,
       height = 5, width = 7.5, dpi = 500)



# Supplemental Figure SXX: PCA of environmental variables -----------------



# Convert environmental data to a matrix
eaa_envData_matrix <- eaa_environmental_data %>%
  select(location_abbr, all_of(eaa_environmental_vars$variable), -starts_with("PC")) %>%
  as.data.frame() %>%
  column_to_rownames("location_abbr") %>%
  as.matrix()

# Run PCA
eaa_envData_pca <- prcomp(x = eaa_envData_matrix)

# Plot
eaa_envData_pcs_df <- eaa_envData_pca$x %>%
  as.data.frame() %>%
  rownames_to_column("location_abbr")

g_envData_pca <- eaa_envData_pcs_df %>%
  ggplot(aes(x = PC1, y = PC2, color = location_abbr)) +
  geom_point() +
  geom_text_repel(aes(label = location_abbr))





# Supplemental Figure S03: dissection of Fst outliers ---------------------


loci_highlight <- fst_outliers %>%
  arrange(desc(Fst)) %>%
  pull(marker)


# Subset this from the list of plots above
allele_geo_list <- sigmar_summary %>%
  filter(marker %in% loci_highlight) %>%
  split(.$marker) %>%
  map("plot2") %>%
  map(1) %>%
  map(~.x + theme_nothing(font_size = base_font_size) + theme(legend.position = .x$theme$legend.position))


## Bioclim variables versus marker alleles
bioclim_allele_list <- sigmar_summary %>%
  filter(marker %in% loci_highlight, !is.na(variable)) %>%
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
      ggplot(aes(x = genotypes_overall, y = biovar, group = genotypes_overall)) +
      geom_jitter(aes(fill = genotypes_overall), pch = 21, size = 1,
                  position = position_jitterdodge(jitter.width = 0.10, dodge.width = 0.75)) +
      geom_point(data = genotype_means, position = position_dodge(0.75), pch = "+", size = 5, color = neyhart_palette()[3]) +
      facet_wrap(~ variable, ncol = 3, scales = "free",
                 labeller = labeller(variable = f_rename_variables)) + # dir = "v" means fill by column
      scale_fill_manual(values = genotype_colors, guide = FALSE) +
      scale_y_continuous(name = NULL, breaks = pretty) +
      scale_x_discrete(name = NULL) +
      theme_grey(base_size = base_font_size) +
      theme(strip.background = element_blank(), panel.grid.minor.x = element_blank(),
            panel.grid.major.x = element_blank(), legend.position = c(0.80, 0.25))

  })


## Combine everything

# List of combined plots
g_combined_list <- full_join(tibble(marker = names(allele_geo_list), p1 = allele_geo_list),
                             tibble(marker = names(bioclim_allele_list), p2 = bioclim_allele_list)) %>%
  mutate_at(vars(p1, p2), ~modify_if(.x = ., .p = sapply(., is.null), ~plot_spacer())) %>%
  mutate(p_combined = map2(p1, p2, ~wrap_plots(.x, .y, ncol = 1, heights = c(1, 0.5))),
         marker = factor(marker, levels = loci_highlight)) %>%
  arrange(marker) %>%
  pull(p_combined)

# Merge
g_combined <- wrap_plots(g_combined_list, ncol = 2)

# Save
ggsave(filename = "figureS03_Fst_outliers.jpg", plot = g_combined, path = fig_dir,
       height = 10, width = 10, dpi = 500)



# Supplemental Figure S04: dissection of SPA outliers ---------------------


loci_highlight <- spa_outliers %>%
  arrange(desc(spa_score)) %>%
  pull(marker)


# Subset this from the list of plots above
allele_geo_list <- sigmar_summary %>%
  filter(marker %in% loci_highlight) %>%
  split(.$marker) %>%
  map("plot2") %>%
  map(1) %>%
  map(~.x + theme_nothing(font_size = base_font_size) + theme(legend.position = .x$theme$legend.position))


## Bioclim variables versus marker alleles
bioclim_allele_list <- sigmar_summary %>%
  filter(marker %in% loci_highlight, !is.na(variable)) %>%
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
      ggplot(aes(x = genotypes_overall, y = biovar, group = genotypes_overall)) +
      geom_jitter(aes(fill = genotypes_overall), pch = 21, size = 1,
                  position = position_jitterdodge(jitter.width = 0.10, dodge.width = 0.75)) +
      geom_point(data = genotype_means, position = position_dodge(0.75), pch = "+", size = 5, color = neyhart_palette()[3]) +
      facet_wrap(~ variable, ncol = 3, scales = "free",
                 labeller = labeller(variable = f_rename_variables)) + # dir = "v" means fill by column
      scale_fill_manual(values = genotype_colors, guide = FALSE) +
      scale_y_continuous(name = NULL, breaks = pretty) +
      scale_x_discrete(name = NULL) +
      theme_grey(base_size = base_font_size) +
      theme(strip.background = element_blank(), panel.grid.minor.x = element_blank(),
            panel.grid.major.x = element_blank(), legend.position = c(0.80, 0.25))

  })


## Combine everything

# List of combined plots
g_combined_list <- full_join(tibble(marker = names(allele_geo_list), p1 = allele_geo_list),
                             tibble(marker = names(bioclim_allele_list), p2 = bioclim_allele_list)) %>%
  mutate_at(vars(p1, p2), ~modify_if(.x = ., .p = sapply(., is.null), ~plot_spacer())) %>%
  mutate(p_combined = map2(p1, p2, ~wrap_plots(.x, .y, ncol = 1, heights = c(1, 0.5))),
         marker = factor(marker, levels = loci_highlight)) %>%
  arrange(marker) %>%
  pull(p_combined)

# Merge
g_combined <- wrap_plots(g_combined_list, ncol = 2)

# Save
ggsave(filename = "figureS04_SPA_outliers.jpg", plot = g_combined, path = fig_dir,
       height = 10, width = 10, dpi = 500)





# Supplemental Table S01: Name and description of environmental variables ---------------

eaa_environmental_vars %>%
  arrange(class, variable) %>%
  rename_all(~str_replace_all(., "_", " ") %>% str_to_title()) %>%
  write_csv(x = ., file = file.path(fig_dir, "table_S01_env_variable_description.csv"))







# Supplemental Table S02: environmental variable PCA -----------------------

pca_env_max_loadings %>%
  mutate_if(is.double, format_numbers) %>%
  rename_all(~str_replace_all(., "_", " ") %>% str_to_title()) %>%
  write_csv(x = ., file = file.path(fig_dir, "table_S02_env_variable_PCA.csv"))



# Supplemental Table S03: all EAA significant markers ---------------------------------


# Get the name and distance of the nearest gene annotation for each marker
all_sigmar_nearby_annotation1 <- all_sigmar_nearby_annotation %>%
  mutate(variable = ifelse(variable == "Fst", "F[ST]", variable)) %>%
  group_by(class, variable, marker) %>%
  do({
    row <- .
    if (nrow(row$within_annotations[[1]]) > 0) {
      anno_name <- cran_gene_ranges[as.numeric(row.names(subset(row$within_annotations[[1]], type == "gene"))),]$Name
      alias <- cran_gene_ranges[as.numeric(row.names(subset(row$within_annotations[[1]], type == "gene"))),]$Alias[[1]]
      within_gene_anno <- tibble(minimum_distance = 0, annotation_name = anno_name,
                                 annotation_alias = setdiff(alias, anno_name))
      start_use <- subset(row$within_annotations[[1]], type == "gene", start, drop = TRUE)

    } else {
      within_gene_anno <- tibble(minimum_distance = as.numeric(NULL), annotation_name = as.character(NULL),
                                 annotation_alias = as.character(NULL))
      start_use <- -1

    }

    # If no annotations, return empty tibble
    if (nrow(row$nearby_annotation[[1]]) == 0) {
      tibble(nNearbyGenes = 0, closest_gene_annotation = as.character(NA),
             alias = as.character(NA), distance = as.numeric(NA))

    } else {

      all_anno <- row$nearby_annotation[[1]] %>%
        filter(type == "gene", start != start_use) %>%
        unnest(cols = c(attributes)) %>%
        mutate(alias = map2(Name, Alias, ~setdiff(str_split(.y, ",")[[1]], .x)),
               alias = map_chr(alias, ~ifelse(length(.x) == 0, as.character(NA), .x))) %>%
        select(minimum_distance = min_distance, annotation_name = Name, annotation_alias = alias)

      bind_rows(within_gene_anno, all_anno) %>%
        mutate(nNearbyGenes = nrow(.)) %>%
        top_n(x = ., n = 1, wt = -minimum_distance) %>%
        select(nNearbyGenes, closest_gene_annotation = annotation_name, alias = annotation_alias, distance = minimum_distance)

    }
  }) %>% ungroup()

# Format the allele frequencies for merging
all_marker_freq1 <- all_marker_freq %>%
  select(marker, category, ref_allele_freq) %>%
  # Convert to MAF
  split(.$marker) %>%
  map_df(~{
    minor <- subset(.x, category == "Wild", ref_allele_freq, drop = TRUE) < 0.5
    if (minor) mutate(.x, maf = ref_allele_freq) else mutate(.x, maf = 1 - ref_allele_freq)
  }) %>%
  select(-ref_allele_freq) %>%
  spread(category, maf)

# add variable and class information
all_marker_freq2 <- all_sigmar %>%
  left_join(., select(eaa_environmental_vars, variable, full_name, class)) %>%
  left_join(., all_marker_freq1) %>%
  mutate(class = ifelse(is.na(class), variable, class)) %>%
  gather(population, maf, Breeding, NativeSelection, Wild) %>%
  select(class, variable, marker, score, population, maf)

# Get the average and range of maf in each population across different association classes
all_marker_freq2 %>%
  group_by(class, population) %>%
  summarize_at(vars(maf), list(mean = mean, min = min, max = max)) %>%
  ungroup() %>%
  as.data.frame()

# class      population      mean         min       max
# 1                F[ST]        Breeding 0.1410613 0.009615385 0.2500000
# 2                F[ST] NativeSelection 0.1008395 0.026666667 0.1933333
# 3                F[ST]            Wild 0.1081081 0.090090090 0.1171171
# 4            geography        Breeding 0.2674247 0.024038462 0.5462963
# 5            geography NativeSelection 0.2622575 0.086419753 0.4444444
# 6            geography            Wild 0.2110682 0.103603604 0.4099099
# 7        precipitation        Breeding 0.2332906 0.000000000 0.4861111
# 8        precipitation NativeSelection 0.2955457 0.020000000 0.4814815
# 9        precipitation            Wild 0.2807207 0.103603604 0.4594595
# 10 principal_component        Breeding 0.2329569 0.028846154 0.6388889
# 11 principal_component NativeSelection 0.2649030 0.086419753 0.6234568
# 12 principal_component            Wild 0.2265122 0.103603604 0.4504505
# 13                soil        Breeding 0.2216200 0.000000000 0.5092593
# 14                soil NativeSelection 0.2203304 0.026666667 0.4814815
# 15                soil            Wild 0.1883943 0.090090090 0.4909910
# 16                 SPA        Breeding 0.3695869 0.057692308 0.7115385
# 17                 SPA NativeSelection 0.3839012 0.246666667 0.6666667
# 18                 SPA            Wild 0.4414414 0.265765766 0.4954955
# 19         temperature        Breeding 0.3469598 0.046296296 0.8221154
# 20         temperature NativeSelection 0.2573489 0.074074074 0.4691358
# 21         temperature            Wild 0.2235657 0.094594595 0.4639640


# Markers with the greatest difference between wild and breeding
all_marker_freq2 %>%
  spread(population, maf) %>%
  arrange(Breeding - Wild) %>%
  as.data.frame()

# These are markers where the MAF in the Wild population is greater
# than that in the Breeding population

# class         variable       marker     score    Breeding NativeSelection       Wild
# 1         precipitation            bio15 S11_32311630 4.3902643 0.043859649      0.26881720 0.38288288
# 2                   SPA              SPA  S10_6176309 2.1116651 0.140350877      0.27419355 0.45045045
# 3                  soil nitrogen_subsoil S09_21132326 4.2430569 0.081818182      0.22413793 0.29729730
# 4                  soil      ocd_subsoil S09_21132326 3.9288414 0.081818182      0.22413793 0.29729730
# 5                   SPA              SPA S06_24012352 2.6171128 0.267543860      0.34408602 0.48198198
# 6                  soil     clay_subsoil S10_28130450 4.6151503 0.000000000      0.13218391 0.21171171
# 7         precipitation            bio14  S02_9985929 3.8055971 0.254385965      0.31720430 0.40090090
# 8         precipitation            bio14 S07_23156165 3.9702409 0.337719298      0.43010753 0.46846847
# 9         precipitation            bio15 S09_17284424 5.4557595 0.008771930      0.10215054 0.13063063
# 10            geography         latitude S02_22935920 4.2805699 0.022727273      0.14367816 0.14414414
# 11        precipitation            bio15   S12_807322 4.1189793 0.000000000      0.01724138 0.11261261
# 12                 soil nitrogen_subsoil S07_13657548 6.2586738 0.004545455      0.14367816 0.11261261
# 13                 soil      ocd_subsoil S07_13657548 4.9050698 0.004545455      0.14367816 0.11261261
# 14                 soil      soc_subsoil S07_13657548 6.3690926 0.004545455      0.14367816 0.11261261
# 15          temperature             bio2 S06_28401995 4.4305959 0.293859649      0.33333333 0.39639640
# 16                 soil     bdod_topsoil S11_15016526 4.4725765 0.000000000      0.06896552 0.09909910
# 17  principal_component              PC3  S01_1383069 4.2005294 0.048245614      0.11827957 0.14414414
# 18                 soil     silt_topsoil  S01_1383069 3.7476058 0.048245614      0.11827957 0.14414414
# 19                F[ST]            F[ST] S02_12681394 0.6596766 0.031818182      0.08620690 0.12162162
# 20                F[ST]            F[ST]  S12_9950985 0.8273571 0.009090909      0.03448276 0.09459459
# 21                 soil     bdod_subsoil  S12_9950985 6.1794434 0.009090909      0.03448276 0.09459459
# 22                 soil     cfvo_subsoil  S12_9950985 4.9311919 0.009090909      0.03448276 0.09459459
# 23                 soil     cfvo_topsoil  S12_9950985 6.2354847 0.009090909      0.03448276 0.09459459
# 24            geography         latitude  S02_5334410 3.8386654 0.346491228      0.45161290 0.40990991
# 25                 soil     silt_topsoil S03_30022546 3.5232702 0.144736842      0.13978495 0.20270270
# 26                 soil nitrogen_subsoil S12_11311413 3.7621700 0.145454545      0.21264368 0.20270270
# 27          temperature             bio4 S04_12490711 4.0928256 0.043859649      0.07526882 0.09459459
# 28                 soil     bdod_subsoil S07_13180166 3.6019107 0.166666667      0.19354839 0.21171171
# 29                 soil    phh2o_subsoil S07_13180166 5.2865427 0.166666667      0.19354839 0.21171171
# 30        precipitation            bio15 S06_32454326 3.7315631 0.144736842      0.14516129 0.18468468
# 31                 soil nitrogen_subsoil  S03_6654650 3.7162006 0.077272727      0.23563218 0.10810811
# 32                 soil      ocd_subsoil  S03_6654650 3.7432647 0.077272727      0.23563218 0.10810811
# 33                 soil     sand_topsoil  S03_6654650 4.6859957 0.077272727      0.23563218 0.10810811
# 34                 soil     silt_subsoil  S03_6654650 5.1657828 0.077272727      0.23563218 0.10810811
# 35                 soil     silt_topsoil  S03_6654650 5.4633502 0.077272727      0.23563218 0.10810811
# 36        precipitation            bio14 S09_18359750 4.0628551 0.090909091      0.10344828 0.12162162
# 37            geography         latitude  S12_9381253 4.3612437 0.136363636      0.22988506 0.16666667
# 38        precipitation            bio15 S09_25633368 4.6274152 0.118181818      0.11494253 0.14414414

# How many markers have elevated MAF in breeding compared to Wild?
all_marker_freq2 %>%
  spread(population, maf) %>%
  filter(! class %in% c("F[ST]", "SPA")) %>%
  distinct(marker, Breeding, NativeSelection, Wild) %>%
  summarize(mafBreed_gt_mafWild = sum(Breeding > Wild),
            mafNS_gt_mafWild = sum(NativeSelection > Wild))

# mafBreed_gt_mafWild mafNS_gt_mafWild
#                 30               38


# What is the genomewide average?
all_marker_freq1 %>%
  summarize(mafBreed_gt_mafWild = sum(Breeding > Wild),
            mafNS_gt_mafWild = sum(NativeSelection > Wild)) %>%
  mutate_all(~. / nrow(snp_info))



all_sigmar %>%
  left_join(., select(all_sigmar_nearby_annotation1, -class)) %>%
  # Add variable classes
  left_join(., select(eaa_environmental_vars, variable, full_name, class)) %>%
  # Add frequency of reference allele in different breeding classes
  left_join(., all_marker_freq1) %>%
  # Format numbers
  mutate_at(vars(Wild, Breeding, NativeSelection), format_numbers) %>%
  # Convert Fst and SPA to their own class; also abbreviate the EAA variables
  mutate(class = ifelse(variable %in% c("F[ST]", "SPA"), variable, paste0("EAA-", str_to_title(str_sub(class, 1, 4)))),
         full_name = ifelse(is.na(full_name), variable, full_name),
         alias2 = closest_gene_annotation,
         closest_gene_annotation = ifelse(is.na(alias), closest_gene_annotation, alias),
         alias = ifelse(is.na(alias), alias, alias2)) %>%
  select(class, variable = full_name, marker, chrom, pos, score, maf_wild = Wild, allele_freq_native = NativeSelection,
         allele_freq_breeding = Breeding, No._nearby_genes = nNearbyGenes, closest_gene_annotation,
         closest_gene_distance = distance, arabidopsis_homolog = alias) %>%
  mutate(class = fct_inorder(class), variable = fct_relevel(variable, "F[ST]", "SPA", after = Inf)) %>%
  arrange(class, variable, marker, chrom, pos) %>%
  rename_all(~str_replace_all(., "_", " ") %>% str_to_title() %>% str_replace(., "Maf", "MAF")) %>%
  # Save as a CSV
  write_csv(x = ., file = file.path(fig_dir, "table_S03_all_significant_markers.csv"), na = "")





