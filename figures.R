# Germplasm collection environmental association
#
# Figures and visualization
#
#

# Load packages
library(raster)
library(tidyverse)
library(readxl)
library(broom)
library(neyhart)
library(patchwork)
library(cowplot)

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




# Read in data ------------------------------------------------------------

load(file.path(data_dir, "germplasm_metadata.RData"))
load(file.path(data_dir, "gc_marker_data.RData"))
load(file.path(data_dir, "germplasm_origin_bioclim_data.RData"))
load(file.path(result_dir, "population_genetics_stats.RData"))
# Load the eGWAS results
load(file.path(result_dir, "eaa_gwas_results.RData"))
load(file.path(result_dir, "egwas_sigmar_analysis.RData"))


# Combine germplasm information
germplasm_info_combined <- germ_meta %>%
  select(accession = formatted_name, variety_designation, category, latitude = origin_latitude,
         longitude = origin_longitude, origin = origin_name) %>%
  mutate(individual = accession)




# Prepare data ------------------------------------------------------------


# Center and scale
eaa_environmental_data1 <- eaa_environmental_data %>%
  select(location_abbr, all_of(eaa_environmental_vars$variable)) %>%
  # mutate_at(vars(contains("bio")), scale, scale = FALSE) %>%
  # mutate_at(vars(contains("bio")), as.numeric) %>%
  mutate(elevation = ifelse(is.na(elevation), 0, elevation))


# Prepare the bioclim data as if it were trait data
germplasm_bioclim_data <- germ_meta %>%
  left_join(., select(eaa_environmental_data, origin_name, location_abbr)) %>%
  # Get the state/province name from the location
  mutate(state = str_sub(string = location_abbr, start = -2, end = -1)) %>%
  filter(category == "wild", !is.na(origin_latitude)) %>%
  select(individual, latitude = origin_latitude, longitude = origin_longitude, state) %>%
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

# Data.frame for native selections and cultivars
geno_hmp_native_cultivars_wild <- rbind(geno_mat_native, geno_mat_cultivar, geno_mat_wild) %>%
  {. + 1} %>% # Convert to ref allele count
  t() %>%
  as.data.frame() %>%
  rownames_to_column("marker") %>%
  inner_join(snp_info, .) %>%
  separate(alleles, c("ref_allele", "alt_allele"), sep = "/")


# Frequency of reference alleles for all markers in all population classes
all_marker_freq <- geno_hmp_native_cultivars_wild %>%
  gather(individual, ref_allele_count, -marker:-alt_allele) %>%
  mutate(class = case_when(
    individual %in% row.names(geno_mat_wild1) ~ "wild",
    individual %in% row.names(geno_mat_native) ~ "native",
    individual %in% row.names(geno_mat_cultivar) ~ "cultivar"
  )) %>%
  group_by(marker, ref_allele, class) %>%
  summarize(ref_allele_freq = mean(ref_allele_count, na.rm = TRUE) / 2, .groups = "drop")






# Marker metadata
marker_metadata <- geno_hmp_wild2 %>%
  distinct_at(vars(marker, chrom, pos, contains("allele")))


# Allele colors
allele_colors <- setNames(neyhart_palette("umn2")[3:4], c("ref_allele", "alt_allele"))

# Subset the native selections
# geno_mat_native <- geno_mat_all[native_sel_meta$germplasm_collection_name,]
geno_mat_native <- geno_mat_all[intersect(row.names(geno_mat_all), native_selection_germplasm_names),]
geno_mat_cultivar <- geno_mat_all[intersect(row.names(geno_mat_all), cultivar_germplasm_names),]


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

# Create a tibble of model parameters to check
egwas_models_df <- tribble(
  ~model, ~K, ~PCs, ~covariates,
  "model1", TRUE, 0, character(),
  "model2", TRUE, 3, character(),
  "model3", TRUE, 0, c("latitude"),
  "model4", TRUE, 3, c("latitude"),
  "model5", TRUE, 0, c("latitude", "longitude"),
  "model6", TRUE, 3, c("latitude", "longitude"),
) %>% mutate(output = list(NULL))


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


  filename <- paste0("egwas_manhattan_", vari, ".jpg")
  ggsave(filename = filename, plot = g_man, path = fig_dir, width = 8, height = 4, dpi = 300)


}


# Sort the manhattan plots
eaa_gwas_manhattan_plot_df <- eaa_gwas_manhattan_plot_df %>%
  mutate(variable = factor(variable, levels = eaa_environmental_vars$variable)) %>%
  arrange(variable) %>%
  mutate(variable = as.character(variable))





# Iterate over markers and create two plots:
# 1. alleles versus bioclim variation
# 2. alleles across geography

# Set seed for reproducible jitter
set.seed(206)
#
sigmar_summary <- egwas_sigmar %>%
  bind_rows(., select(filter(gwas_out_best_model, marker == "S07_30582060", variable == "bio1"), variable, best_model = model, marker, p_value)) %>%
  group_by(variable, marker) %>%
  do({
    row <- .
    marker_i <- row$marker
    var_i <- row$variable
    var_i_fullname <- subset(eaa_environmental_vars, variable == var_i, full_name, drop = TRUE)


    # Get the mlmm output for this variable
    biovar1 <- select(germplasm_bioclim_data, individual, biovar1 = all_of(var_i))

    # Subset the bioclim data for this variable
    bioclim_data_i <- germplasm_bioclim_data %>%
      select(individual, lat = latitude, long = longitude, location_abbr) %>%
      left_join(., select(eaa_environmental_data, location_abbr, biovar = all_of(var_i)), by = "location_abbr") %>%
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
      left_join(., select(eaa_environmental_data, location_abbr, biovar = all_of(var_i)), by = "location_abbr") %>%
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
      spread(allele, freq)



    # Create genotype colors
    genotype_colors <- colorRampPalette(colors = allele_colors)(3) %>%
      setNames(., c(paste0(ref_allele, ref_allele), paste0(ref_allele, alt_allele), paste0(alt_allele, alt_allele)))


    # Plot marker allele versus bioclim
    g1 <- ggplot(data = bioclim_data_i, aes(x = genotypes_overall, y = biovar, group = genotype)) +
      # geom_boxplot(width = 0.5, alpha = 0.5, outlier.shape = NA) +
      geom_violin(aes(fill = genotype), alpha = 0.5) +
      geom_jitter(width = 0.1, size = 0.5) +
      xlab(marker_i) +
      ylab(str_wrap(var_i_fullname, width = 30)) +
      scale_fill_manual(values = genotype_colors, guide = FALSE) +
      neyhart::theme_genetics()

    g1_alt <- ggplot(data = bioclim_data_i, aes(x = genotypes_overall, y = biovar1, group = genotype)) +
      geom_boxplot(width = 0.5, alpha = 0.5, outlier.shape = NA) +
      geom_jitter(width = 0.1, size = 0.5) +
      xlab(marker_i) +
      ylab(str_wrap(var_i_fullname, width = 30)) +
      neyhart::theme_genetics()

    colors <- structure(RColorBrewer::brewer.pal(n = 11, name = "RdBu"), class = "palette")[c(3, 6, 9)]
    colors <- RColorBrewer::brewer.pal(n = 9, name = "Blues")[c(2, 5, 8)]

    # Plot frequency of alleles in different geographies
    # Plot marker alleles across geographies
    # g2 <- g_map +
    #   geom_point(data = bioclim_data_i, aes(x = long, y = lat, fill = genotype), shape = 21,
    #              position = position_jitter(0.1, 0.1)) +
    #   scale_fill_manual(values = colors, name = "Marker allele") +
    #   theme(legend.position = c(0.8, 0.2), legend.background = element_rect(color = "black"),
    #         legend.margin = margin(3,3,3,3))


    # Rename allele colors
    allele_colors1 <- allele_colors
    names(allele_colors1) <- c(ref_allele, alt_allele)

    g2 <- g_map +
      scatterpie::geom_scatterpie(data = bioclim_data_i2, aes(x = long, y = lat, group = location_abbr),
                                  cols = unlist(subset(marker_metadata, marker == marker_i, c(ref_allele, alt_allele))),
                                  size = 0.5, pie_scale = 1, alpha = 0.75, legend_name = "Allele") +
      scale_fill_manual(values = allele_colors1, name = "Allele\nfrequencies") +
      # coord_map(projection = "bonne", lat0 = mean(g_map$coordinates$limits$y),
      #           xlim = g_map$coordinates$limits$x, ylim = g_map$coordinates$limits$y)
      coord_fixed(xlim = g_map$coordinates$limits$x, ylim = g_map$coordinates$limits$y)

    # Alternative map with individual points and intermediate colors for hets
    g2_alt <- g_map +
      geom_jitter(data = bioclim_data_i, aes(color = genotype), width = 0.25, height = 0.25, size = 1) +
      # scale_color_manual(values = genotype_colors, name = "Marker\ngenotype") +
      # theme(legend.position = c(0, 0), legend.justification = c(0, 0),
      #       legend.box.background = element_rect(fill = "white", colour = NA), legend.box.margin = margin(2, 5, 2, 5))
      scale_color_manual(values = genotype_colors, name = "Marker genotype",
                         guide = guide_legend(override.aes = list(size = 2))) +
      theme(legend.position = "bottom", legend.justification = "left")


    # Calculate the frequency of each allele in the wild, native, and cultivar populations
    ref_allele_freq <- tribble(
      ~population, ~allele_type, ~allele, ~frequency,
      "wild", "ref", ref_allele, mean(geno_mat_wild1[,marker_i] + 1) / 2,
      "native_selections", "ref", ref_allele,  mean(geno_mat_native[,marker_i] + 1) / 2,
      # "cultivars", "ref", ref_allele, mean(geno_mat_cultivar[,marker_i] + 1) / 2
    )
    allele_freq <- ref_allele_freq %>%
      bind_rows(., mutate(ref_allele_freq, allele_type = "alt", allele = alt_allele, frequency = 1 - frequency)) %>%
      mutate(population = fct_inorder(population),
             allele = factor(allele, levels = c("A", "C", "G", "T")))

    # Add an annotation for Fst
    fst_out <- subset(marker_fst_wild_native, marker == marker_i, Fstp, drop = TRUE)
    # Quantile of this value
    fst_out_quant <- mean(marker_fst_wild_native$Fstp <= fst_out)

    # Plot this too
    g3 <- ggplot(data = allele_freq, aes(x = population, y = frequency, fill = allele)) +
      geom_col() +
      geom_text(aes(label = round(frequency, 3)), position = position_stack(vjust = 0.5)) +
      labs(subtitle = paste0("Fst: ", format_numbers(x = fst_out, signif.digits = 3), " (",
                             format_numbers(x = fst_out_quant * 100, signif.digits = 3), "%)")) +
      theme_genetics()


    # Return a tibble
    tibble(plot1 = list(g1),
           plot1_alt = list(g1_alt),
           plot2 = list(g2),
           plot2_alt = list(g2_alt),
           plot3 = list(g3),
           allele_freq = list(allele_freq)
    )


  }) %>% ungroup() %>%
  left_join(., select(eaa_environmental_vars, -class))



# Create a bunch of plots
pb <- progress::progress_bar$new(total = nrow(sigmar_summary))



for (i in seq_len(nrow(sigmar_summary))) {

  # Determine the left and right plots
  left_plot <- sigmar_summary$plot2_alt[[i]] +
    theme(legend.position = "none")
  right_plot <- sigmar_summary$plot1[[i]]

  combined_plot1 <- plot_grid(left_plot, right_plot, rel_widths = c(1, 0.4))


  # Save
  filename <- paste0("egwas_summary_", sigmar_summary$variable[i], "-marker-", sigmar_summary$marker[i], "-",
                     sigmar_summary$class[i], ".jpg")
  # Set seed for reproducible jitter
  set.seed(206)
  ggsave(filename = filename, plot = combined_plot1, path = fig_dir, width = 6, height = 2.5, dpi = 300)

  pb$tick()

}












# Figure 1: geographic origin, population structure, ld decay -------------

# Plot origin of wild germplasm

wild_germ_meta1 <- germ_meta %>%
  rename(latitude = origin_latitude, longitude = origin_longitude) %>%
  filter(!is.na(latitude)) %>%
  rename(long = longitude, lat = latitude)

# Summarize number of genotyped by site
wild_germ_meta_origin_summ <- wild_germ_meta1 %>%
  group_by(long, lat) %>%
  summarize(nSamples = n(), .groups = "drop")

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
# PCA of both wild and native selections
# K_pca <- prcomp(x = K_all)

# Gather the eigenvector
pc_eigenvec <- K_pca$rotation %>%
  as.data.frame() %>%
  rownames_to_column("individual") %>%
  as_tibble() %>%
  left_join(., bind_rows(tibble(category = "wild", individual = wild_germplasm_names), tibble(category = "native_selections", individual = native_selection_germplasm_names)))

# Get the eigenvalues
eigenvals <- K_pca$sdev^2

# Calculate proportion of explained variance
pc_varexp <- eigenvals / sum(eigenvals)
names(pc_varexp) <- paste0("PC", seq_along(eigenvals))

# Cumulative sum
plot(cumsum(pc_varexp))
head(pc_varexp)
sum(pc_varexp[1:3])


# Create a plotting df
pca_plot_df <- pc_eigenvec %>%
  left_join(., distinct(germ_meta, individual, origin_name)) %>%
  inner_join(., select(eaa_environmental_data, origin_name, location_abbr, state))

# Calculate the median position of each color and plot a label for that
group_label <- pca_plot_df %>%
  group_by(state) %>%
  summarize_at(vars(contains("PC")), median)

# Manual range
x_limits <- c(-0.15, 0.25)
y_limits <- c(-0.2, 0.45)

# add Origin and plot
pc_plot <- pca_plot_df %>%
  ggplot(aes(x = PC1, y = PC2, color = state)) +
  geom_point(aes(alpha = location_abbr, linetype = individual)) +
  ggrepel::geom_label_repel(data = group_label, aes(label = state), fill = NA, key_glyph = "point", size = 2,
                            max.overlaps = 15, point.padding = unit(10, "line"), box.padding = unit(0.2, "line")) +
  scale_color_discrete(guide = FALSE) +
  scale_x_continuous(name = paste0("PC1 (", round(pc_varexp["PC1"] * 100, 3), "%)"), breaks = pretty, limits = x_limits,
                     labels = format_numbers) +
  scale_y_continuous(name = paste0("PC2 (", round(pc_varexp["PC2"] * 100, 3), "%)"), breaks = pretty, limits = y_limits,
                     labels = format_numbers) +
  scale_alpha_discrete(guide = FALSE) +
  theme_genetics() +
  theme(legend.position = c(0.99, 0.99), legend.justification = c(0.99, 0.99))

# Save this
ggsave(filename = "snp_pca_wild_germplasm.jpg", plot = pc_plot, path = fig_dir, width = 5, height = 4, dpi = 1000)


# Plot PCA with the state color in a separate legend
pc_plot1 <- pca_plot_df %>%
  ggplot(aes(x = PC1, y = PC2, color = state)) +
  geom_point(aes(alpha = location_abbr, linetype = individual)) +
  scale_color_discrete(name = "Origin state", guide = guide_legend(ncol = 2)) +
  scale_x_continuous(name = paste0("PC1 (", round(pc_varexp["PC1"] * 100, 3), "%)"), breaks = pretty, labels = format_numbers) +
  scale_y_continuous(name = paste0("PC2 (", round(pc_varexp["PC2"] * 100, 3), "%)"), breaks = pretty, labels = format_numbers) +
  scale_alpha_discrete(guide = FALSE) +
  theme_genetics() +
  theme(legend.position = c(0.5, 0.5), legend.justification = c(0, 0))


# Save this
ggsave(filename = "snp_pca_wild_germplasm1.jpg", plot = pc_plot1, path = fig_dir, width = 5, height = 4, dpi = 1000)





###

# Modify the fitted LD values
ld_wild_decay1 <- unnest(ld_wild_decay, predictions)
ld_wild_df1 <- ld_wild_df %>%
  # rename(d = bp_dist) %>%
  filter(d <= max(ld_wild_decay1$d))

## Plot LD decay

# Plot LD decay without filtering distance

g_ld_decay_all_chrom1 <- ld_wild_df1 %>%
  ggplot(aes(x = d / 1000, y = r2)) +
  geom_point(size = 0.3, color = "gray") +
  geom_line(data = ld_wild_decay1, aes(x = d / 1000, y = yhat), lwd = 1) +
  scale_x_continuous(name = "Distance (kbp)", breaks = pretty) +
  scale_y_continuous(name = expression("Linkage disequilibrium ("*r^2*")")) +
  facet_wrap(~ LG) +
  theme_genetics()


# Save
ggsave(filename = "ld_decay_all_chrom_1000kbp.jpg", plot = g_ld_decay_all_chrom1, path = fig_dir,
       width = 5, height = 4, dpi = 1000)


# Plot LD decay with filtering distance

# Modify the fitted LD values
ld_wild_decay2 <- unnest(ld_wild_decay, predictions) %>%
  filter(d <= 2e5)
ld_wild_df2 <- ld_wild_df %>%
  # rename(d = bp_dist) %>%
  filter(d <= 2e5)


g_ld_decay_all_chrom2 <- ld_wild_df2 %>%
  ggplot(aes(x = d / 1000, y = r2)) +
  geom_point(size = 0.3, color = "gray") +
  geom_line(data = ld_wild_decay2, aes(x = d / 1000, y = yhat), lwd = 1) +
  scale_x_continuous(name = "Distance (kbp)", breaks = pretty) +
  scale_y_continuous(name = expression("Linkage disequilibrium ("*r^2*")")) +
  facet_wrap(~ LG) +
  theme_genetics()



# Save
ggsave(filename = "ld_decay_all_chrom_200kbp.jpg", plot = g_ld_decay_all_chrom2, path = fig_dir,
       width = 5, height = 4, dpi = 1000)


## Plot chromosome 1 low distance versus full

g_ld_decay_chrom1_long <- ld_wild_df1 %>%
  filter(LG == "01") %>%
  ggplot(aes(x = d / 1000, y = r2)) +
  geom_point(size = 0.3, color = "gray") +
  geom_line(data = subset(ld_wild_decay1, LG == "01"), aes(x = d / 1000, y = yhat), lwd = 1) +
  scale_x_continuous(name = "Chromosome 1 marker distance (kbp)", breaks = pretty, expand = expansion(mult = c(0.01, 0.05))) +
  scale_y_continuous(name = expression("Linkage disequilibrium ("*r^2*")"), expand = expansion(mult = c(0.01, 0.05))) +
  theme_genetics()

g_ld_decay_chrom1_short <- ld_wild_df2 %>%
  filter(LG == "01") %>%
  ggplot(aes(x = d / 1000, y = r2)) +
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

# Save
ggsave(filename = "ld_decay_combined_chrom1.jpg", plot = combined_ld_plot, path = fig_dir,
       width = 5, height = 4, dpi = 1000)


## Plot whole genome SPA and Fst
##
##

library(slider)

# Combine SPA with Fst
spa_fst <- full_join(spa_out, marker_fst_wild)

# Correlation
cor(spa_fst$spa_score, spa_fst$Fst)
plot(spa_fst$spa_score, spa_fst$Fst)



spa_fst_slide <- spa_fst %>%
  # Calculate sliding window of values
  mutate_at(vars(spa_score, Fst), list(slide = ~slide_dbl(., mean, .before = 4, .after = 4, .step = 5))) %>%
  filter(!is.na(spa_score_slide))


spa_fst %>%
  ggplot(aes(x = pos / 1e6)) +
  geom_point(aes(y = spa_score), color = "red") +
  geom_point(aes(y = -Fst), color = "blue") +
  facet_grid(~ chrom, space = "free_x", scales = "free_x", switch = "x") +
  scale_x_continuous(name = NULL, breaks = aggregate(pos ~ chrom, spa_fst, median)$pos) +
  theme_genetics(8) +
  theme(strip.placement = "outside", strip.background.x = element_blank(),
        panel.spacing = unit(0, "line"))

spa_fst_slide %>%
  ggplot(aes(x = pos / 1e6)) +
  # geom_point(aes(y = spa_score_slide), color = "red") +
  # geom_point(aes(y = -Fst_slide), color = "blue") +
  geom_line(aes(y = spa_score_slide), color = "red") +
  geom_line(aes(y = -Fst_slide), color = "blue") +
  facet_grid(~ chrom, space = "free_x", scales = "free_x", switch = "x") +
  scale_x_continuous(name = NULL, breaks = aggregate(pos ~ chrom, spa_fst, median)$pos) +
  theme_genetics(8) +
  theme(strip.placement = "outside", strip.background.x = element_blank(),
        panel.spacing = unit(0, "line"))



# Any overlap?
intersect(marker_fst_wild_outlier$marker, spa_outliers$marker)








# Combine pop structure with LD
g_popstr_ld <- plot_grid(pc_plot, combined_ld_plot, ncol = 1, labels = subfigure_labels[1:2],
                         label_size = 9, align = "h", axis = "tblr")

# Save
ggsave(filename = "figure2_popstr_chrom1_ld.jpg", plot = g_popstr_ld, path = fig_dir,
       width = 3.5, height = 6, dpi = 1000)


# ## Add pairwise genetic distance data?
#
#
# # Convert the pairwise genetic distance df to a matrix
#
# pairwise_dist <- pairwise_dist_data %>%
#   select(contains("individual"), distance) %>%
#   spread(individual1, distance) %>%
#   as.data.frame() %>%
#   column_to_rownames("individual2") %>%
#   as.dist()
#
# # Pairwise distance of individuals
# pairwise_geo_dist <- germplasm_bioclim_data %>%
#   select(individual, latitude, longitude) %>%
#   crossing(., ., .name_repair = tidyr_legacy) %>%
#   mutate(gcd = geosphere::distGeo(p1 = select(., longitude, latitude), p2 = select(., longitude1, latitude1))) %>%
#   select(individual, individual1, gcd) %>%
#   # Filter
#   filter_at(vars(contains("individual")), all_vars(. %in% unique(pairwise_dist_data$individual1))) %>%
#   spread(individual, gcd) %>%
#   as.data.frame() %>%
#   column_to_rownames("individual1") %>%
#   as.dist()
#
# # Conduct mantel test
# mantel_out <- vegan::mantel(xdis = pairwise_dist, ydis = pairwise_geo_dist)
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






# Figure 3: Temperature related associations -----------------------------


## First show the locus on Chromosome 4 associated with multiple worldclim temperature
## variables

loci_highlight <- c("S04_24963382", "S07_30582060")

locus1_summary <- egwas_sigmar %>%
  filter(marker %in% loci_highlight, str_detect(variable, "bio")) %>%
  # add environmental variable full names
  left_join(., select(eaa_environmental_vars, -class))

# Create a highlight window for this snp
locus1_highlight_rect <- snp_info %>%
  filter(marker %in% loci_highlight) %>%
  mutate(left = pos - 1000, right = pos + 1000)


# Create a merged manhattan plot
#
# Subset manhattan plots for biovars associated with this locus
manhattan_list <- subset(eaa_gwas_manhattan_plot_df, variable %in% locus1_summary$variable, plot, drop = TRUE) %>%
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
g_manhattan_combine <- cowplot::plot_grid(plotlist = manhattan_list, ncol = 1,
                                          rel_heights = c(1, rep(0.7, length(manhattan_list)-2), 0.85))
g_manhattan_combine1 <- cowplot::plot_grid(y_axis_label, g_manhattan_combine, rel_widths = c(0.03, 1))

# Save
ggsave(filename = "temperature_ROI_manhattan_merge.jpg", plot = g_manhattan_combine1, path = fig_dir,
       height = 5, width = 3.5, dpi = 1000)




## Plot the geographic distribution of the alleles at these SNPs
sigmar_allele_geo_dist <- sigmar_summary %>%
  select(-variable, -full_name) %>%
  split(.$marker) %>%
  map_df(head, 1)

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
  filter(marker %in% loci_highlight, str_detect(variable, "bio")) %>%
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
sigmar_summary_temperature_list <- map(names(allele_geo_list), ~{
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
  plot_grid(top_plot1, bottom_plot, ncol = 1, rel_heights = c(1, 0.8))

})

names(sigmar_summary_temperature_list) <- names(allele_geo_list)

# Save plots
for (i in seq_along(sigmar_summary_temperature_list)) {
  marker <- names(sigmar_summary_temperature_list)[i]
  filename <- paste0("example_associations_marker_", marker, ".jpg")
  ggsave(filename = filename, plot = sigmar_summary_temperature_list[[i]], path = fig_dir,
         height = 3, width = 4, dpi = 1000)

}



### Combine with the manhattan plots


# Merge using patchwork
combined_plot <- g_manhattan_combine1 + (sigmar_summary_temperature_list[[1]] / sigmar_summary_temperature_list[[2]]) +
  plot_layout(widths = c(0.65, 1)) +
  plot_annotation(tag_levels = "A") +
  theme(plot.tag = element_text(size = base_font_size + 2))

# Save
ggsave(filename = "figure3_temperature_associations.jpg", plot = combined_plot, path = fig_dir,
       height = 6, width = 9, dpi = 1000)





# Display the frequency of each allele in each population
allele_geo_list$S04_24963382$layers[[4]]$data %>%
  group_by(location_abbr, genotype) %>%
  count() %>%
  ungroup() %>%
  spread(genotype, n) %>%
  mutate_at(vars(-location_abbr), ~ifelse(is.na(.), 0, .)) %>%
  mutate(minor_count = (select(., 2)[[1]] * 2) + (select(., 3)[[1]]),
         major_count = (select(., 4)[[1]] * 2) + (select(., 3)[[1]]),
         minor_freq = minor_count / (minor_count + major_count))

# Display the frequency of each allele in each population
allele_geo_list[[2]]$layers[[4]]$data %>%
  group_by(location_abbr, genotype) %>%
  count() %>%
  ungroup() %>%
  spread(genotype, n) %>%
  mutate_at(vars(-location_abbr), ~ifelse(is.na(.), 0, .)) %>%
  mutate(minor_count = (select(., 2)[[1]] * 2) + (select(., 3)[[1]]),
         major_count = (select(., 4)[[1]] * 2) + (select(., 3)[[1]]),
         minor_freq = minor_count / (minor_count + major_count)) %>%
  arrange(minor_freq)



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


# Frequency of alleles
subset(all_marker_freq, marker %in% loci_highlight)



# Nearby genes for these loci?
egwas_sigmar_nearby_annotation1 %>%
  filter(marker %in% loci_highlight) %>%
  distinct(marker, nearby_annotation, closest_gene_dist) %>%
  unnest() %>%
  filter(type == "gene") %>%
  filter(end_distance == min(end_distance)) %>%
  unnest(attributes) %>%
  as.data.frame()






# Figure 4: associations related to soil variables --------------------------




## First show the locus on Chromosome 4 associated with multiple worldclim temperature
## variables

loci_highlight <- c("S02_38027635")

locus1_summary <- egwas_sigmar %>%
  filter(marker %in% loci_highlight) %>%
  # add environmental variable full names
  left_join(., select(eaa_environmental_vars, -class))

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

allele_freq <- sigmar_allele_geo_dist %>%
  filter(marker %in% loci_highlight) %>%
  select(marker, allele_freq)



## Bioclim variables versus marker alleles
bioclim_allele_list <- sigmar_summary %>%
  filter(marker %in% loci_highlight, variable %in% locus1_summary1$variable) %>%
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
combined_plot <- g_manhattan_combine1 + sigmar_summary_soil_list +
  plot_layout(widths = c(0.65, 1)) +
  plot_annotation(tag_levels = "A") +
  theme(plot.tag = element_text(size = base_font_size + 2))

# Save
ggsave2(filename = "figure4_soil_associations", plot = sigmar_summary_soil_list$S02_38027635,
       height = 3, width = 3.5)


bioclim_allele_list$S02_38027635$data %>%
  group_by(variable, genotype) %>%
  summarize(biovar_mean = mean(biovar), nInd = n(), .groups = "drop") %>%
  filter(map_lgl(str_split(genotype, ""), ~n_distinct(.) == 1)) %>%
  mutate(genotype_class = ifelse(nInd == min(nInd), "minor", "major")) %>%
  select(-genotype, -nInd) %>%
  spread(genotype_class, biovar_mean) %>%
  mutate(diff = minor - major)

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
  distinct(marker, nearby_annotation, closest_gene_dist) %>%
  unnest() %>%
  filter(type == "gene") %>%
  filter(end_distance == min(end_distance)) %>%
  unnest(attributes) %>%
  as.data.frame()






# Table S1: Name and description of environmental variables ---------------

eaa_environmental_vars %>%
  filter(str_detect(variable, "prec|tmax|tmin|tmean", negate = TRUE)) %>%
  filter(! variable %in% c("elevation", "latitude", "longitude", "PC1", "PC2", "PC3")) %>%
  arrange(class, variable) %>%
  rename_all(~str_replace_all(., "_", " ") %>% str_to_title()) %>%
  write_csv(x = ., file = file.path(fig_dir, "table_S01_env_variable_description.csv"))




# Table S02: All significant associations ---------------------------------

# Add the full names of the environmental variables
egwas_sigmar %>%
  filter(str_detect(variable, "prec|tmax|tmin|tmean", negate = TRUE)) %>%
  select(-model, -class, -contains("score")) %>%
  left_join(., eaa_environmental_vars) %>%
  select(variable, variable_full_name = full_name, names(.)) %>%
  # Save as a CSV
  write_csv(x = ., file = file.path(fig_dir, "table_S02_egwas_significant_associations.csv"))





















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
