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
# Read in the GRIN germplasm metadata
grin_germ_meta <- read_excel(path = file.path(cran_dir, "Breeding/Germplasm/GRINGermplasm/Vm_GRIN_search_20210329.xlsx")) %>%
  select(-starts_with("..."), -IMAGE) %>%
  rename_all(~str_replace_all(tolower(.), " ", "_")) %>%
  # Parse lat/long
  separate(coordinates, c("latitude", "longitude"), sep = ", ") %>%
  mutate_if(is.character, parse_guess) %>%
  # Sanity check for coordinates
  mutate(longitude = ifelse(sign(longitude) == 1, -1 * longitude, longitude)) %>%
  mutate_at(vars(source_type, improvement_status), tolower) %>%
  mutate(accession = str_replace_all(accession, " ", "-"),
         source = "grin-ncgr")

# Combine germplasm information
germplasm_info_combined <- germ_meta %>%
  select(accession = formatted_name, variety_designation, category, latitude = origin_latitude,
         longitude = origin_longitude, origin = origin_name) %>%
  mutate(individual = accession) %>%
  bind_rows(., select(grin_germ_meta, accession, individual = name, origin,
                      category = improvement_status, latitude, longitude, availability))

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



# Plot germplasm origin - this is a base plot for later plots ---------------------------------------------------

## Plot sampling location

# Create a box for x and y axes
xlim <- range(pretty(germplasm_bioclim_data$long))
ylim <- range(pretty(germplasm_bioclim_data$lat))

# Get the map data for canada
canada <- rnaturalearth::ne_states(country = "canada") %>%
  broom::tidy(x = ., region = "name_en") %>%
  mutate(group = as.numeric(as.factor(group))) %>%
  mutate(area = "canada")

# Load US state data
us_states <- map_data("state") %>%
  # Adjust the groups in the states
  mutate(group = group + max(canada$group),
         area = "usa_state")

north_america_mapdata <- bind_rows(us_states, canada)

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




# Create manhattan plots for each association -----------------------------




# Subset GWAS results for the best model
gwas_out_best_model <- gwas_out %>%
  inner_join(., subset(p_lambda, best_model == "*", c(model, variable), drop = TRUE))


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

  scores1 <- gwas_out %>%
    filter(variable == vari, model == subset(p_lambda, variable == vari & best_model == "*", model, drop = TRUE)) %>%
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
  ggsave(filename = filename, plot = g_man, path = fig_dir, width = 8, height = 4, dpi = 1000)


}


# Sort the manhattan plots
eaa_gwas_manhattan_plot_df <- eaa_gwas_manhattan_plot_df %>%
  mutate(variable = factor(variable, levels = eaa_environmental_vars$variable)) %>%
  arrange(variable) %>%
  mutate(variable = as.character(variable))





# Compare significant marker alleles with bioclim variables ---------------

# Extract the significant markers from the multi-locus mixed model output
egwas_sigmar <- eaa_gwas_best_model_signif_marker %>%
  unnest(sigmar) %>%
  mutate(class = "egwas_sigmar")

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

# In addition, find the markers with the three most extreme p-values
outlier_sigmar <- gwas_out %>%
  inner_join(., select(filter(p_lambda, best_model == "*"), variable, model)) %>%
  split(.$variable) %>%
  map_df(~arrange(., p_value) %>% head(3)) %>%
  select(marker, variable) %>%
  mutate(class = "outlier") %>%
  # Remove those markers found in the mlmm
  anti_join(x = ., y = distinct(egwas_sigmar, marker, variable))

all_sigmar <- bind_rows(egwas_sigmar, outlier_sigmar)


# Tranpose the genotype hmp object
geno_hmp_wild2 <- geno_hmp_wild1 %>%
  separate(col = alleles, c("ref_allele", "alt_allele"), sep = "/") %>%
  mutate(marker1 = paste0(marker, "_", ref_allele, "_", alt_allele)) %>%
  filter(marker %in% unique(union(egwas_sigmar$marker, outlier_sigmar$marker))) %>%
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



# Iterate over markers and create two plots:
# 1. alleles versus bioclim variation
# 2. alleles across geography

# Set seed for reproducible jitter
set.seed(206)
#
sigmar_summary <- all_sigmar %>%
  group_by(variable, marker, class) %>%
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
      geom_jitter(data = bioclim_data_i, aes(color = genotype), width = 0.25, height = 0.25) +
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
  left_plot <- sigmar_summary$plot2_alt[[i]]
  right_plot <- sigmar_summary$plot1[[i]]

  ## Combine plots
  layout <- c(
    area(t = 1, l = 1, b = 10, r = 8),
    area(t = 6, l = 7, b = 10, r = 9)
  )

  combined_plot1 <- plot_grid(left_plot) +
    plot_grid(right_plot)+
    plot_layout(design = layout)


  # Save
  filename <- paste0("egwas_summary_", sigmar_summary$variable[i], "-marker-", sigmar_summary$marker[i], "-",
                     sigmar_summary$class[i], ".jpg")
  ggsave(filename = filename, plot = combined_plot1, path = fig_dir, width = 7, height = 4, dpi = 500)

  pb$tick()

}







# Compare eGWAS results with Sambada and SPA --------------------------------------

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











# Highlight example results for presentations -----------------------------

biovar_select <- c("bio1", "phh2o_topsoil", "sand_topsoil")


## First show the locus on Chromosome 4 associated with multiple worldclim temperature
## variables

locus1_highlight <- "S04_24963382"

locus1_summary <- sigmar_summary %>%
  filter(term == locus1_highlight, str_detect(variable, "bio")) %>%
  # Only keep bios 1, 2, 6, 7
  filter(variable %in% paste0("bio", c(1,2,6,7))) %>%
  arrange(parse_number(variable))

# Create a highlight window for this snp
locus1_highlight_rect <- snp_info %>%
  filter(marker == locus1_highlight) %>%
  mutate(left = pos - 10000, right = pos + 10000)

# Create a merged manhattan plot
#
# Subset manhattan plots for biovars associated with this locus
manhattan_list <- subset(eaa_gwas_manhattan_plot_df, variable %in% locus1_summary$variable, plot, drop = TRUE) %>%
  map(., ~. + geom_vline(data = locus1_highlight_rect, aes(xintercept = pos), lwd = 5, color = "grey75", alpha = 0.2) +
        theme(axis.title.y = element_blank())) %>%
  map(~{.x$layers[[2]]$aes_params$size <- 1.5; .x}) %>% # Modify point size
  modify_at(.x = ., .at = -length(.), ~. + theme(axis.ticks.x = element_blank(), strip.text.x = element_blank())) %>%
  modify_at(.x = ., .at = 1, ~ . + theme(legend.direction = "horizontal")) %>%
  modify_at(.x = ., .at = -1, ~. + theme(legend.position = "none"))

# Text grob
y_axis_label <- grid::textGrob(label = expression(-log[10](p-value)), rot = 90,
                               gp = grid::gpar(fontsize = 16))

# Combine these
g_manhattan_combine <- cowplot::plot_grid(plotlist = manhattan_list, ncol = 1,
                                          rel_heights = c(1, rep(0.7, length(manhattan_list)-2), 0.85))
g_manhattan_combine1 <- cowplot::plot_grid(y_axis_label, g_manhattan_combine, rel_widths = c(0.03, 1))

# Save
ggsave(filename = paste0(locus1_highlight, "_manhattan_merge.jpg"), plot = g_manhattan_combine1, path = fig_dir,
       height = 8, width = 10, dpi = 1000)



## Merge the snp alleles vs bioclim data plot
effect_plot_list <- locus1_summary$plot1 %>%
  map(., ~. + theme_genetics(base_size = 16) + theme(axis.title.x = element_blank())) %>%
  map(., ~. + ylab(str_wrap(str_trim(str_split(string = .x$labels$y, pattern = "\\(")[[1]][1]), 20))) %>%
  map(~{
    .x$layers[[2]]$aes_params$size <- 1
    .x
  })

g_plot1_combine <- cowplot::plot_grid(plotlist = effect_plot_list, nrow = 1, align = "hv")

ggsave(filename = paste0(locus1_highlight, "_biovar_effect_plot.jpg"), plot = g_plot1_combine, path = fig_dir,
       height = 3, width = 12, dpi = 1000)


# Save a larger version of the allele frequency distribution plot
g_map_use <- locus1_summary$plot2[[1]] +
  theme_void(base_size = 20) +
  # labs(subtitle = "Spatial allele frequencies") +
  theme(legend.position = c(0.85, 0.1), legend.justification = c(1,0))

ggsave(filename = paste0(locus1_highlight, "_spatial_frequency.jpg"), plot = g_map_use, path = fig_dir,
       height = 4, width = 12, dpi = 1000)



##
## Second look at soil pH and sand content associations
##
loci_summary <- sigmar_summary %>%
  filter(variable %in% biovar_select,
         str_detect(variable, "bio", negate = TRUE),
         class == "mlmm") %>%
  left_join(., mlmm_out_list[biovar_select] %>% map(1) %>% imap_dfr(~mutate(.x, variable = .y))) %>%
  group_by(variable) %>%
  top_n(x = ., n = 1, wt = -p_value) %>%
  # # Pull the fdr thresholds
  # left_join(., mutate(eaa_gwas_manhattan_plot_df, fdr = map_dbl(plot, ~min(.x$layers[[1]]$data$fdr_score, na.rm = TRUE)))) %>%
  # # Highlight regions above -log10(p) >= 4
  # filter(-log10(p_value) >= fdr) %>%
  # Manual filter
  filter(
    (variable == "sand_topsoil" & str_detect(term, "S02|S03|S06")) | variable == "phh2o_topsoil",
    p_value < 0.0025
  )


# Subset the mahattan plots
#
#
# Create highlights of the loci
loci2_highlight_rect <- snp_info %>%
  inner_join(., select(loci_summary, variable, term), by = c("marker" = "term")) %>%
  mutate(left = pos - 10000, right = pos + 10000) %>%
  group_by(variable) %>%
  nest() %>%
  ungroup()

manhattan_list <- eaa_gwas_manhattan_plot_df %>%
  inner_join(., loci2_highlight_rect) %>%
  mutate(plot = map2(.x = plot, .y = data, ~{
    .x1 <- .x + geom_vline(data = .y, aes(xintercept = pos), lwd = 5, color = "grey75", alpha = 0.2) +
      theme(axis.title.y = element_blank(), legend.direction = "horizontal")
    .x1$labels$subtitle <- str_replace(string = .x1$labels$subtitle, pattern = "Ph", replacement = "pH")
    .x1$layers[[2]]$aes_params$size <- 1.5 # Modify point size
    .x1
  }))

# Save these
for (i in seq_len(nrow(manhattan_list))) {
  filename <- paste0(manhattan_list$variable[i], "_manhattan_toplot.jpg")
  ggsave(filename = filename, plot = manhattan_list$plot[[i]], path = fig_dir,
         width = 10, height = 3, dpi = 1000)
}


# Subset the effect plots
effect_plot_list <- setNames(object = loci_summary$plot1, nm = loci_summary$variable) %>%
  map(., ~. + theme_genetics(base_size = 16) + theme(axis.title.x = element_blank())) %>%
  map(., ~. + ylab(str_wrap(str_trim(str_split(string = .x$labels$y, pattern = "\\(")[[1]][1]), 20))) %>%
  map(., ~{
    .x$labels$y <- str_replace(string = .x$labels$y, pattern = "Ph", replacement = "pH")
    .x$layers[[2]]$aes_params$size <- 1
    .x
  })

# Subset the map plots
map_plot_list <- setNames(object = loci_summary$plot2, nm = loci_summary$variable) %>%
  map(., ~{
    .x + theme_void(base_size = 16) +
      # labs(subtitle = "Spatial allele frequencies") +
      theme(legend.position = c(0.90, 0.05), legend.justification = c(1,0))
  })

# Merge plots
loci_effect_map_merge <- map2(.x = map_plot_list, .y = effect_plot_list, ~{
  cowplot::plot_grid(.x, .y, rel_widths = c(1, 0.6))
})

# Save these
for (i in seq_along(loci_effect_map_merge)) {
  filename <- paste0(names(loci_effect_map_merge)[i], "_map_effect_toplot.jpg")
  ggsave(filename = filename, plot = loci_effect_map_merge[[i]], path = fig_dir,
         width = 10, height = 2, dpi = 1000)
}




# Other summary statistics
all_assoc <- filter(sigmar_summary, class == "mlmm")
nrow(all_assoc)
n_distinct(all_assoc$variable)
n_distinct(all_assoc$term)

xtabs(~ chrom, snp_info, subset = marker %in% all_assoc$term)




























