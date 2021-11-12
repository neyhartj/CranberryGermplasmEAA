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



# Read in data ------------------------------------------------------------

load(file.path(data_dir, "germplasm_metadata.RData"))
load(file.path(data_dir, "gc_marker_data.RData"))
load(file.path(data_dir, "germplasm_origin_bioclim_data.RData"))
load(file.path(result_dir, "population_genetics_stats.RData"))



# Read in the germplasm metadata
germ_meta <- mutate(germ_meta, source = "rutgers")

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
         longitude = origin_longitude, origin = origin_name, source) %>%
  mutate(individual = accession) %>%
  bind_rows(., select(grin_germ_meta, accession, individual = name, origin,
                      category = improvement_status, latitude, longitude, source, availability))

# Read in world clim data
# Read in already downloaded worldclim data
worldclim_dat <- getData(name = "worldclim", download = FALSE, var = "bio", res = "5",
                         path = file.path(env_dir, "WorldClim"))






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
ggsave(filename = "figure1_cranberry_wild_germplasm_origins1.tiff", plot = g_map1, path = fig_dir,
       width = 8, height = 4, dpi = 1000)



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
  scale_x_continuous(name = paste0("PC1 (", round(pc_varexp["PC1"] * 100, 3), "%)"), breaks = pretty, limits = x_limits) +
  scale_y_continuous(name = paste0("PC2 (", round(pc_varexp["PC2"] * 100, 3), "%)"), breaks = pretty, limits = y_limits) +
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
  scale_x_continuous(name = paste0("PC1 (", round(pc_varexp["PC1"] * 100, 3), "%)"), breaks = pretty) +
  scale_y_continuous(name = paste0("PC2 (", round(pc_varexp["PC2"] * 100, 3), "%)"), breaks = pretty) +
  scale_alpha_discrete(guide = FALSE) +
  theme_genetics() +
  theme(legend.position = c(0.5, 0.5), legend.justification = c(0, 0))


# Save this
ggsave(filename = "snp_pca_wild_germplasm1.jpg", plot = pc_plot1, path = fig_dir, width = 5, height = 4, dpi = 1000)


# Modify the fitted LD values
ld_wild_decay1 <- unnest(ld_wild_decay, predictions)
ld_wild_df1 <- ld_wild_df %>%
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
  theme(axis.title = element_blank(), plot.background = element_rect(fill = alpha("white", 0), color = "black"))

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

# Combine pop structure with LD
g_popstr_ld <- plot_grid(pc_plot, combined_ld_plot, ncol = 1, labels = subfigure_labels[1:2],
                         label_size = 9)
# Save
ggsave(filename = "figure2_popstr_chrom1_ld.tiff", plot = g_popstr_ld, path = fig_dir,
       width = 3.5, height = 6, dpi = 1000)








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



