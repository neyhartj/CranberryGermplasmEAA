# Germplasm collection environmental association
#
# Tests for loci with large changes in allele frequency between wild cranberry and the native selections
#

# Load packages
library(sommer)
library(tidyverse)
library(snps)
library(hierfstat)
library(LDheatmap)
library(snpStats)
library(ape)


# Load the startup script
source("startup.R")


# Read in the marker data
load(file.path(data_dir, "gc_marker_data.RData"))
# Read in the germplasm metadata
load(file.path(data_dir, "germplasm_metadata.RData"))
load(file.path(data_dir, "germplasm_origin_bioclim_data.RData"))

germ_meta1 <- germ_meta %>%
  left_join(., distinct(eaa_environmental_data, origin_name, location_abbr))

# Read in variety origin information
native_sel_meta <- read_csv(file = file.path(data_dir, "native_selection_metadata.csv")) %>%
  rename(state = original_planting_state)

# Subset the geno mat for wild and natives
geno_mat_wild_native <- geno_mat_unfiltered[,colnames(geno_mat_wild)]



# Calculate LD decay ------------------------------------------------------

snp_info1 <- snp_info %>%
  filter(marker %in% colnames(geno_mat_wild)) %>%
  select(Locus = marker, LG = chrom, Position = pos) %>%
  as.data.frame()


# Calculate LD decay for each population
ld_wild <- LD.decay(markers = geno_mat_wild, map = snp_info1)
# Combine
ld_wild_df <- ld_wild$all.LG %>%
  rownames_to_column("id") %>%
  mutate(id = str_split(id, "\\."), LG = map_chr(id, 1)) %>%
  as_tibble() %>%
  filter(LG %in% unique(snp_info1$LG))


# plot(r2 ~ d, data = ld_wild_df, subset = LG == "01")

# Plot
g_ld_wild <- ld_wild_df %>%
  filter(LG == "08") %>%
  ggplot(aes(x = d, y = r2)) +
  geom_point(size = 0.5) +
  facet_wrap(~ LG) +
  theme_classic()



## Fit the non-linear models of Abecasis et al. 2001 as coded by Anderson et al 2016:

#   Next, we have to set the starting guesses for our parameters to minimize
#   the RMSE. These can't be exactly 0, since that would produce infinity.
#   nls() expects them to be in a list
start <- list(dlow=0.05, dhigh=0.95, gens=5000)

#   Set the lower bounds. These should be in the same order as the starting
#   values.
lower <- c(0.001, 0.001, 1)
upper <- c(0.999, 0.999, 100000)

# Use a recombination rate of 1 cM / 1 Mb
recomb_rate <- 1
# Scaling factor
scaling <- (1e6 * 100 * 10)

# Combine dfs and calculate LD decay for both natives and wild
ld_wild_decay <- ld_wild_df %>%
  group_by(LG) %>%
  do({
    df <- .
    ld <- df$r2
    d <- df$d
    # Calculate recombination fraction (Mb * (cM/Mb)) / 100 cM;
    # Use this same formula in the model for r
    r <- (((d / 1e6) * recomb_rate) / 100)

    # Fit the model
    ld_decay_model <- nls(formula = ld ~ dlow + (dhigh - dlow) * ((1 - (d / scaling))^gens) * 1,
                          start = start, control = nls.control(maxiter = 100, warnOnly = T),
                          lower = lower, upper = upper, algorithm = "port")

    # Predict
    xvals <- seq(0, 1e6, 500)
    yhat <- predict(object = ld_decay_model, list(d = xvals))
    predictions <- tibble(d = xvals, d_kbp = xvals / 1000,  yhat = yhat)

    model_summary <- summary(ld_decay_model)
    dlow <- model_summary$parameters["dlow", "Estimate"]
    dhigh <- model_summary$parameters["dhigh", "Estimate"]
    gens <- model_summary$parameters["gens", "Estimate"]

    # Calculate the half-life
    half <- yhat[1]/2
    half_scaled <- 1 - 10^(log10((half - dlow)/(dhigh - dlow))/gens)
    half_bp <- (half_scaled * scaling) / recomb_rate

    # Calculate a moving average of LD by d
    # First cut d into windows of size 1 kbp
    window_size <- 10000
    d_range <- range(pretty(d))
    d_seq <- seq(0, max(d_range), by = window_size)
    d_cut <- cut(x = d, breaks = d_seq, right = FALSE)
    LD_moving_avg <- map_dbl(split(x = ld, f = d_cut), mean)
    # Tibble
    LD_moving_df <- tibble(d = head(d_seq, -1), LD_avg = LD_moving_avg)

    # Return a tibble of the parameters
    tibble(predictions = list(predictions), half_bp, LD_moving = list(LD_moving_df))

  }) %>% ungroup()

# Mean and range half-life
ld_wild_decay %>%
  summarize(mean_half_life = mean(half_bp), min_half_life = min(half_bp), max_half_life = max(half_bp))

# When does LD fall below 0.3, 0.1?
ld_wild_decay %>%
  select(LG, predictions) %>%
  crossing(., ld_threshold = c(0.3, 0.1)) %>%
  unnest(predictions) %>%
  filter(yhat < ld_threshold) %>%
  group_by(LG, ld_threshold) %>%
  top_n(x = ., n = 1, wt = -d) %>%
  group_by(ld_threshold) %>%
  summarize(mean_d = mean(d), min_d = min(d), max_d = max(d))



# Calculate pairwise distance between individuals -------------------------

pairwise_dist <- dist.gene(x = geno_mat_wild, method = "percentage")

# Convert to df
pairwise_dist_df <- broom::tidy(pairwise_dist) %>%
  rename(individual1 = item1, individual2 = item2)

# Add location information
pairwise_dist_df1 <- pairwise_dist_df %>%
  left_join(., distinct(germ_meta1, individual, location_abbr), by = c("individual1" = "individual")) %>%
  rename(location1 = location_abbr) %>%
  left_join(., distinct(germ_meta1, individual, location_abbr), by = c("individual2" = "individual")) %>%
  rename(location2 = location_abbr)

# Summarize within-location pairwise distances
pairwise_dist_df1 %>%
  filter(location1 == location2) %>%
  group_by(location1) %>%
  summarize(mean_distance = mean(distance), min_dist = min(distance), max_distance = max(distance)) %>%
  left_join(., aggregate(individual ~ location_abbr, data = germ_meta1, FUN = n_distinct), by = c("location1" = "location_abbr")) %>%
  arrange(desc(mean_distance))

# Rename
pairwise_dist_data <- pairwise_dist_df1



# Run SPA -----------------------------------------------------------------

# Run the run_SPA.sh script from R
#
# Breakdown:
# "bash -c" is used to invoke a bash command
# Then pass the command that one would normally run on the command line: "sh script.sh"
#
system('bash -c "sh R/popGenStats/run_spa.sh /"')

# Read in the spa output
spa_out <- read_tsv(file.path(result_dir, "wild_cranberry_knownGeo_spa.out"),
                    col_names = c("chrom", "marker", "cM", "pos", "ref_allele", "alt_allele", "slope_a", "slope_b",
                                  "slope_c", "spa_score"))



# Calculate Fst between wild and native -----------------------------------

# Create a data.frame to store the results
marker_fst_wild_native <- geno_hmp_wild %>%
  select(marker, chrom, pos) %>%
  mutate(fst_wild_native = as.numeric(NA))

# Create a df specifying population
pop_df <- data.frame(individual = row.names(geno_mat_unfiltered), stringsAsFactors = FALSE) %>%
  mutate(pop = ifelse(individual %in% row.names(geno_mat_wild), 0, 1)) %>%
  mutate(pop = ifelse(pop == 0, "wild", "native")) %>%
  column_to_rownames("individual")

# Format the data at once
geno_mat_wild_native_hierfFormat <- apply(X = geno_mat_wild_native + 1, MARGIN = 2, FUN = function(snp) {
  as.numeric(str_replace_all(snp, c("2" = "22", "1" = "12", "0" = "11")))
})
row.names(geno_mat_wild_native_hierfFormat) <- row.names(geno_mat_wild_native)

# Combine pop_df and geno_mat_wild_native_hierfFormat
data_for_fstat <- cbind(pop_df, as.data.frame(geno_mat_wild_native_hierfFormat))

# Calculate basic stats
basic_stats_out <- basic.stats(data = data_for_fstat)


# Create a summary data.frame per locus
marker_fst_wild_native <- geno_hmp_wild %>%
  select(marker, chrom, pos) %>%
  left_join(., basic_stats_out$perloc %>% as.data.frame() %>% rownames_to_column("marker"))

# What is the 99% percentile of Fst values?
fst_cutoff <- quantile(x = marker_fst_wild_native$Fstp, probs = 0.999)

# Plot
marker_fst_wild_native %>%
  mutate(chrom = parse_number(chrom), even_chrom = chrom %% 2 == 0) %>%
  ggplot(aes(x = pos / 1e6, y = Fstp, color = even_chrom)) +
  geom_hline(yintercept = fst_cutoff, lty = 2) +
  geom_point() +
  scale_color_discrete(guide = FALSE) +
  facet_grid(~ chrom, scales = "free_x", switch = "x") +
  theme_classic() +
  theme(panel.spacing.x = unit(0, "line"))


# List SNPs that exceed the empirical threshold
marker_fst_wild_native %>%
  filter(Fstp >= fst_cutoff) %>%
  arrange(desc(Fstp))



# Save everything ---------------------------------------------------------



# Save the results
save("marker_fst_wild_native", "spa_out", "ld_wild_decay", "ld_wild_df", "pairwise_dist_data",
     file = file.path(result_dir, "population_genetics_stats.RData"))














