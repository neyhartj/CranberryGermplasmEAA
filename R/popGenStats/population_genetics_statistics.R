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


# Split by chromosome and create a list of SnpMatrix objects
ld_wild <- snp_info1 %>%
  split(.$LG) %>%
  map(~{
    # Subset markers from this chrom; create a snpMatrix object
    markers <- .x$Locus
    snpMat <- as(geno_mat_wild[,markers, drop = FALSE] + 1, "SnpMatrix")
    # Physical distances
    bp_dist <- .x$Position

    # Calculate LD
    ld_D <- LDheatmap(gdat = snpMat, genetic.distances = bp_dist, distances = "physical", LDmeasure = "D",
                      add.map = FALSE, flip = TRUE)
    # Calculate LD
    ld_r <- LDheatmap(gdat = snpMat, genetic.distances = bp_dist, distances = "physical", LDmeasure = "r",
                      add.map = FALSE, flip = TRUE)

    # Return a tidy version of the distance matrix
    ld_df1 <- t(ld_D$LDmatrix) %>%
      as.dist() %>%
      broom::tidy() %>%
      rename(marker1 = item1, marker2 = item2, D = distance) %>%
      left_join(., select(snp_info1, Locus, Position), by = c("marker1" = "Locus")) %>%
      left_join(., select(snp_info1, Locus, Position), by = c("marker2" = "Locus")) %>%
      mutate(bp_dist = abs(Position.x - Position.y)) %>%
      select(marker1, marker2, bp_dist, D)

    # Return a tidy version of the distance matrix
    t(ld_r$LDmatrix) %>%
      as.dist() %>%
      broom::tidy() %>%
      rename(marker1 = item1, marker2 = item2, r2 = distance) %>%
      left_join(ld_df1, ., by = c("marker1", "marker2"))

  })


# Combine
ld_wild_df <- imap_dfr(ld_wild, ~mutate(.x, LG = .y))


# Plot
g_ld_wild <- ld_wild_df %>%
  filter(LG == "08") %>%
  ggplot(aes(x = bp_dist, y = r2)) +
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
    # ld <- df$D
    d <- df$bp_dist

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

    # ## plot?
    # plot(ld ~ d)
    # points()

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

ld_wild_decay

# Mean and range half-life
ld_wild_decay %>%
  summarize(mean_half_life = mean(half_bp), min_half_life = min(half_bp), max_half_life = max(half_bp))

# # A tibble: 1 x 3
# mean_half_life min_half_life max_half_life
# 1         21997.         8996.        40233


# When does LD fall below 0.3, 0.2, 0.1?
ld_wild_decay %>%
  select(LG, predictions) %>%
  crossing(., ld_threshold = c(0.3, 0.2, 0.1)) %>%
  unnest(predictions) %>%
  filter(yhat < ld_threshold) %>%
  group_by(LG, ld_threshold) %>%
  top_n(x = ., n = 1, wt = -d) %>%
  group_by(ld_threshold) %>%
  summarize(mean_d = mean(d), min_d = min(d), max_d = max(d))

# ld_threshold mean_d min_d max_d
# 1          0.1  39625 17000 63000
# 2          0.2  16500  7500 23500
# 3          0.3   3875     0  8000



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





# Calculate the number of segregating sites per population ------------------

# Split individuals by population
indiv_per_pop <- germ_meta1 %>%
  filter(individual %in% rownames(geno_mat_wild)) %>%
  split(.$location_abbr) %>%
  map("individual")

pop_seg_sites <- indiv_per_pop %>%
  imap_dfr(~{
    mat1 <- geno_mat_wild[.x,, drop = FALSE]
    ref_allele_freq <- colMeans(mat1 + 1) / 2
    tibble(location_abbr = .y, nSegSites = sum(ref_allele_freq < 1 & ref_allele_freq > 0)) %>%
      mutate(propSegsites = nSegSites / ncol(mat1), popSize = nrow(mat1))
  })




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




# Calculate Fst among all wild individuals --------------------------------

# Create a data.frame to store the results
marker_fst <- geno_hmp_wild %>%
  select(marker, chrom, pos) %>%
  mutate(fst = as.numeric(NA))

# DF with population indicators
pop_df <- germ_meta %>%
  filter(individual %in% row.names(geno_mat_wild)) %>%
  select(individual, population = origin_name)

# Format the data at once
geno_mat_wild_hierfFormat <- apply(X = geno_mat_wild + 1, MARGIN = 2, FUN = function(snp) {
  as.numeric(str_replace_all(snp, c("2" = "22", "1" = "12", "0" = "11")))
})
row.names(geno_mat_wild_hierfFormat) <- row.names(geno_mat_wild)

# Combine pop_df and geno_mat_wild_native_hierfFormat
data_for_fstat <- cbind(pop_df, as.data.frame(geno_mat_wild_hierfFormat)) %>%
  remove_rownames() %>%
  column_to_rownames("individual")

# Calculate basic stats
basic_stats_out <- basic.stats(data = data_for_fstat, diploid = TRUE)

# Use the Weir and Cockeram Fst estimate
wc_fst <- wc(ndat = data_for_fstat, diploid = TRUE)

# Create a summary data.frame per locus
marker_fst_wild <- geno_hmp_wild %>%
  select(marker, chrom, pos) %>%
  # left_join(., basic_stats_out$perloc %>% as.data.frame() %>% rownames_to_column("marker"))
  left_join(., rownames_to_column(wc_fst$per.loc, "marker") %>% rename(Fst = FST, Fis = FIS))

# What is the 99% percentile of Fst values?
fst_cutoff <- quantile(x = marker_fst_wild$Fst, probs = 0.999)


# Plot
g_fst <- marker_fst_wild %>%
  mutate(chrom = parse_number(chrom), even_chrom = chrom %% 2 == 0) %>%
  ggplot(aes(x = pos, y = Fst, color = even_chrom)) +
  geom_hline(yintercept = fst_cutoff, lty = 2) +
  geom_point() +
  scale_color_discrete(guide = FALSE) +
  scale_x_continuous() +
  facet_grid(~ chrom, scales = "free_x", switch = "x") +
  theme_classic() +
  theme(panel.spacing.x = unit(0, "line"))

# Plotly
# plotly::ggplotly(g_fst)


# List SNPs that exceed the empirical threshold
marker_fst_wild_outlier <- marker_fst_wild %>%
  filter(Fst >= fst_cutoff) %>%
  arrange(desc(Fst))







# Save everything ---------------------------------------------------------



# Save the results
save("marker_fst_wild", "marker_fst_wild_outlier", "spa_out", "ld_wild_decay", "ld_wild_df", "pairwise_dist_data",
     file = file.path(result_dir, "population_genetics_stats.RData"))














