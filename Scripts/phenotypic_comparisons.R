# Germplasm environmental association
# 
# Visualize phenotypic distributions per sub-population
# Compare native selections and wild material
# 

# Load packages
library(tidyverse)
library(readxl)

# Set directories
data_dir <- here::here("data")
results_dir <- here::here("results")

# External directories
# CranberryLab directory
path_split <- str_split(string = here::here(), pattern = "/")[[1]]
cran_dir <- paste0(path_split[seq_len(str_which(string = path_split, pattern = "CranberryLab"))], collapse = "/")
# Directory of metadata
meta_dir <- file.path(cran_dir, "Breeding/prior2021/Populations/")



# Read in data ------------------------------------------------------------

# Read in the phenotypic data
pheno_dat <- read_csv(file.path(data_dir, "gc_phenotypic_data.csv"))

# Read in population metadata
metadata <- read_excel(path = file.path(meta_dir, "germplasmCollectionMetadata.xlsx"), sheet = "all_germplasm")

# Create abbreviations for the origin names
origin_names_abbr <- metadata %>%
  filter(!is.na(origin_name)) %>% 
  distinct(origin_name, origin_latitude, origin_longitude) %>%
  mutate(origin_name1 = str_remove(origin_name, "west branch, ")) %>% 
  separate(origin_name1, c("location", "state"), sep = ", ") %>% 
  mutate(origin_abbr = abbreviate(location, 3),
         origin_abbr = paste0(toupper(origin_abbr), state))

# Assign origin and category metadata to the phenotypic data
pheno_dat1 <- left_join(pheno_dat, select(metadata, individual, category, origin_name))


## Number of individuals in each category
pheno_dat1 %>%
  distinct(individual, category) %>%
  group_by(category) %>%
  summarize(n = n())


# Visualize ---------------------------------------------------------------

# Vector of ordered traits
trait_order <- c("fruit_yield", "fruit_weight", "fruit_rot_percent", "brix", "tacy", "titrat_acidity")

# Filter phenotypic data for critical traits
pheno_dat2 <- pheno_dat1 %>%
  filter(trait %in% trait_order)

# Cultivar checks
cultivar_pheno_dat <- pheno_dat2 %>%
  filter(individual %in% c("CRIMSON_QUEEN", "DEMORANVILLE", "MULLICA_QUEEN"))

# Compare phenotypic distributions of native selections and wilds
pheno_dat2 %>%
  filter(category %in% c("wild", "native_selection", "production_selection")) %>%
  ggplot(aes(x = value, fill = category)) + 
  geom_density(alpha = 0.5) +
  geom_vline(data = cultivar_pheno_dat, aes(xintercept = value, color = individual)) +
  facet_wrap(~ trait, scales = "free")


## PCA of traits
pheno_mat_scaled <- pheno_dat2 %>%
  select(individual, trait, value) %>%
  spread(trait, value) %>%
  as.data.frame() %>%
  column_to_rownames("individual") %>%
  as.matrix() %>%
  scale()

# Impute missing with 0
pheno_mat_scaled[is.na(pheno_mat_scaled)] <- 0

pca_out <- prcomp(x = pheno_mat_scaled, retx = TRUE)
# PCA variance explained
pc_varexp <- pca_out$sdev / sum(pca_out$sdev)
names(pc_varexp) <- paste0("PC", seq_along(pc_varexp))

# Retried results; add category designation
pca_out_tidy <- pca_out$x %>% 
  as.data.frame() %>%
  rownames_to_column("individual") %>%
  left_join(., select(metadata, individual, category))

# Filter and plot
pca_out_tidy1 <- pca_out_tidy %>%
  filter(category %in% c("wild", "native_selection", "production_selection"))

g_pheno_pc <- pca_out_tidy1 %>%
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(aes(color = category), size = 1.5) +
  # add cultivars
  geom_point(data = inner_join(pca_out_tidy, distinct(cultivar_pheno_dat, individual)),
             aes(shape = individual), size = 2) +
  scale_x_continuous(name = paste0("PC (", round(pc_varexp, 3)["PC1"] * 100, "%)"), breaks = pretty) +
  scale_y_continuous(name = paste0("PC (", round(pc_varexp, 3)["PC2"] * 100, "%)"), breaks = pretty) +
  labs(subtitle = paste0("Phenotype PCA\n(", paste0(colnames(pheno_mat_scaled), collapse = ", "), ")")) +
  theme_bw()



# Examine phenotypic differences across locations -------------------------

pheno_dat_origin <- pheno_dat2 %>% 
  filter(!is.na(origin_name)) %>%
  left_join(., origin_names_abbr)

# Plot
pheno_dat_origin %>%
  ggplot(aes(x = value, fill = origin_abbr)) + 
  geom_density(alpha = 0.5) +
  # geom_vline(data = cultivar_pheno_dat, aes(xintercept = value, color = individual)) +
  facet_wrap(~ trait, scales = "free")

# Summarize
pheno_dat_origin_summ <- pheno_dat_origin %>%
  group_by(trait, origin_abbr) %>%
  summarize(mean = mean(value, na.rm = TRUE),
            sd = sd(value, na.rm = TRUE),
            n = n()) %>%
  ungroup() %>%
  mutate(se = sd / sqrt(n), lower = mean - se, upper = mean + se,
         origin_abbr1 = paste0(origin_abbr, " (n=", n, ")"))


# Order and plot
g_origin_pheno <- pheno_dat_origin_summ %>%
  split(.$trait) %>%
  map_df(~mutate(.x, order = rank(mean))) %>%
  mutate(order = as.factor(order), trait = factor(trait, levels = trait_order)) %>%
  ggplot(aes(x = order, y = mean, color = origin_abbr1)) +
  geom_errorbar(aes(ymin = lower, ymax = upper)) +
  geom_point() +
  facet_wrap(~ trait, scales = "free") +
  scale_x_discrete(labels = NULL) +
  scale_y_continuous(breaks = pretty) +
  theme_bw()

plotly::ggplotly(g_origin_pheno)


# Alternative plot
g_origin_pheno1 <- pheno_dat_origin %>%
  group_by(trait, location) %>%
  mutate(median = median(value, na.rm = TRUE)) %>%
  ungroup() %>%
  split(.$trait) %>%
  map_df(~mutate(.x, order = rank(median))) %>%
  mutate(order = as.factor(order), trait = factor(trait, levels = trait_order)) %>%
  ggplot(aes(x = order, y = value, fill = origin_abbr, group = origin_abbr)) +
  geom_boxplot(alpha = 0.5) +
  facet_wrap(~ trait, scales = "free") +
  scale_x_discrete(labels = NULL) +
  scale_y_continuous(breaks = pretty) +
  theme_bw()

plotly::ggplotly(g_origin_pheno1)


## Very clear variation of phenotypes among origins