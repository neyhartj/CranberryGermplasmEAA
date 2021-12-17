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


# Read in the base data
load(file.path(data_dir, "population_metadata_and_genotypes.RData"))
load(file.path(data_dir, "germplasm_origin_bioclim_data.RData"))

germ_meta1 <- pop_metadata %>%
  left_join(., select(eaa_environmental_data, location_of_origin, location_abbr))


# Calculate LD decay ------------------------------------------------------

# Read in the vcf- this is PHASED data
vcf <- vcfR::read.vcfR(file = file.path(data_dir, "wild_cranberry_cleaned_genotypes.vcf.gz"))

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
    snpMat <- vcfR2SnpMatrix(obj = vcf[which(vcf@fix[,"ID"] %in% markers),], phased = TRUE)
    # Physical distances
    bp_dist <- .x$Position

    # Calculate LD
    ld_D <- LDheatmap(gdat = snpMat$data, genetic.distances = bp_dist, distances = "physical", LDmeasure = "D",
                      add.map = FALSE, flip = TRUE)
    # Calculate LD
    ld_r <- LDheatmap(gdat = snpMat$data, genetic.distances = bp_dist, distances = "physical", LDmeasure = "r",
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
    ld_df2 <- t(ld_r$LDmatrix) %>%
      as.dist() %>%
      broom::tidy() %>%
      rename(marker1 = item1, marker2 = item2, r2 = distance) %>%
      left_join(ld_df1, ., by = c("marker1", "marker2"))


    # Simple correlation matrix
    ld_cor2 <- (cor(geno_mat_wild[,markers])^2) %>%
      as.dist() %>%
      broom::tidy() %>%
      rename(marker1 = item1, marker2 = item2, cor2 = distance) %>%
      left_join(ld_df2, ., by = c("marker1", "marker2"))

  })


# Combine
ld_wild_df <- imap_dfr(ld_wild, ~mutate(.x, LG = .y))


# Plot
g_ld_wild <- ld_wild_df %>%
  filter(LG == "08") %>%
  ggplot(aes(x = bp_dist, y = cor2)) +
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

# Fit decay models
# Combine dfs and calculate LD decay for both natives and wild
ld_wild_decay <- ld_wild_df %>%
  group_by(LG) %>%
  do({
    df <- .
    # ld <- df$r2
    # ld <- df$D
    ld <- df$cor2
    d <- df$bp_dist

    # Calculate recombination fraction (Mb * (cM/Mb)) / 100 cM;
    # Use this same formula in the model for r
    r <- (((d / 1e6) * recomb_rate) / 100)

    # Fit the model
    ld_decay_model <- nls(formula = ld ~ dlow + (dhigh - dlow) * ( (1 - (d / scaling) * 1)^gens),
                          start = start, control = nls.control(maxiter = 1000, warnOnly = T),
                          lower = lower, upper = upper, algorithm = "port", trace = FALSE)

    # Predict
    xvals <- seq(0, 1e6, 500)
    yhat <- predict(object = ld_decay_model, list(d = xvals))
    predictions <- tibble(d = xvals, d_kbp = xvals / 1000,  yhat = yhat)

    # ## plot?
    # plot(ld ~ d)
    # points()

    model_summary <- summary(ld_decay_model)
    dlow_out <- model_summary$parameters["dlow", "Estimate"]
    dhigh_out <- model_summary$parameters["dhigh", "Estimate"]
    gens_out <- model_summary$parameters["gens", "Estimate"]

    # Calculate the half-life
    half <- yhat[1]/2
    half_scaled <- 1 - 10^(log10((half - dlow_out)/(dhigh_out - dlow_out))/gens_out)
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

# Quick plot
ld_wild_decay %>%
  subset(LG == "01") %>%
  unnest(predictions) %>%
  plot(yhat ~ d, ., ylim = c(0, 1))

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

# location_abbr nSegSites propSegsites popSize
# 1 VNLME               294       0.0601       1
# 2 MTPPA              1509       0.309        1
# 3 LBCME              1517       0.310        2
# 4 FIENY              1678       0.343        1
# 5 PTCNY              2334       0.477        2
# 6 ONCWI              2816       0.576        4
# 7 AMGNY              3004       0.614        3
# 8 SNNMA              3313       0.677        4
# 9 STCNB              3382       0.691        3
# 10 FLNNY              3455       0.706        6
# 11 FIWNY              3641       0.744        4
# 12 LWSDE              3763       0.769        8
# 13 CRLMI              4117       0.842        9
# 14 IBSPNJ             4210       0.861       14
# 15 SSQPA              4280       0.875        6
# 16 CRMWV              4444       0.909        8
# 17 PLLMA              4739       0.969       35


# Run SPA -----------------------------------------------------------------

# Run the run_SPA.sh script from R
#
# Breakdown:
# "bash -c" is used to invoke a bash command
# Then pass the command that one would normally run on the command line: "sh script.sh"
#
system('bash -c "sh Scripts/popGenStats/run_spa.sh /"')

# Read in the spa output
spa_out <- read_tsv(file.path(result_dir, "wild_cranberry_knownGeo_spa.out"),
                    col_names = c("chrom", "marker", "cM", "pos", "ref_allele", "alt_allele", "slope_a", "slope_b",
                                  "slope_c", "spa_score"))




# Calculate Fst among all wild individuals --------------------------------

# DF with population indicators
pop_df <- pop_metadata %>%
  filter(individual %in% row.names(geno_mat_wild)) %>%
  select(individual, population = location_of_origin)

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
# basic_stats_out <- basic.stats(data = data_for_fstat, diploid = TRUE)

# Use the Weir and Cockeram Fst estimate
wc_fst <- wc(ndat = data_for_fstat, diploid = TRUE)

# Create a summary data.frame per locus
global_marker_fst <- geno_hmp_wild %>%
  select(marker, chrom, pos) %>%
  # left_join(., basic_stats_out$perloc %>% as.data.frame() %>% rownames_to_column("marker"))
  left_join(., rownames_to_column(wc_fst$per.loc, "marker") %>% rename(Fst = FST, Fis = FIS))

# Save global Fst
(global_fst <- wc_fst$FST)

## Calculate Fst while removing singleton locations
wc_fst2 <- data_for_fstat %>%
  inner_join(., subset(as.data.frame(xtabs(~ population, data_for_fstat)), Freq > 1)) %>%
  wc(ndat = ., diploid = TRUE)

# Create a summary data.frame per locus
global_marker_fst2 <- geno_hmp_wild %>%
  select(marker, chrom, pos) %>%
  # left_join(., basic_stats_out$perloc %>% as.data.frame() %>% rownames_to_column("marker"))
  left_join(., rownames_to_column(wc_fst2$per.loc, "marker") %>% rename(Fst = FST, Fis = FIS))

(global_fst2 <- wc_fst2$FST)

# Calculate pairwise Fst between populations
pairwise_data_for_fst <- pop_df %>%
  distinct(population) %>%
  crossing(., ., .name_repair = tidyr_legacy) %>%
  filter(population != population1) %>%
  mutate(data = map2(population, population1, ~subset(data_for_fstat, population %in% c(.x, .y))))

pairwise_wc_fst <- pairwise.WCfst(dat = data_for_fstat, diploid = TRUE)

# Neighbor join
pairwise_wc_fst_nj1 <- pairwise_wc_fst_nj <- nj(X = pairwise_wc_fst)
# Adjust tip labels
pairwise_wc_fst_nj1$tip.label <-  abbreviate(str_remove_all(pairwise_wc_fst_nj1$tip.label, ", "))
# Plot
plot(x = pairwise_wc_fst_nj1, type = "u")



# Plot
g_fst <- global_marker_fst2 %>%
  mutate(chrom = parse_number(chrom), even_chrom = chrom %% 2 == 0) %>%
  ggplot(aes(x = pos, y = Fst, color = even_chrom)) +
  geom_hline(yintercept = quantile(global_marker_fst2$Fst, 0.999), lty = 2) +
  geom_point() +
  scale_color_discrete(guide = FALSE) +
  scale_x_continuous() +
  facet_grid(~ chrom, scales = "free_x", switch = "x") +
  theme_classic() +
  theme(panel.spacing.x = unit(0, "line"))

# Plotly
# plotly::ggplotly(g_fst)


# List SNPs that exceed the empirical threshold
global_marker_fst_wild_outlier <- global_marker_fst_wild %>%
  filter(Fst >= fst_cutoff) %>%
  arrange(desc(Fst))



# Pairwise Fst between wild/native/breeding -------------------------------

# DF with population indicators
pop_df <- pop_metadata %>%
  filter(individual %in% row.names(geno_mat_all)) %>%
  select(individual, population = category)

# Format the data at once
geno_mat_all_hierfFormat <- apply(X = geno_mat_all + 1, MARGIN = 2, FUN = function(snp) {
  as.numeric(str_replace_all(snp, c("2" = "22", "1" = "12", "0" = "11")))
})
row.names(geno_mat_all_hierfFormat) <- row.names(geno_mat_all)

# Combine pop_df and geno_mat_wild_native_hierfFormat
data_for_fstat <- cbind(pop_df, as.data.frame(geno_mat_all_hierfFormat)) %>%
  remove_rownames() %>%
  column_to_rownames("individual")

# Calculate basic stats
# basic_stats_out <- basic.stats(data = data_for_fstat, diploid = TRUE)

# Use the Weir and Cockeram Fst estimate
wc_fst_WildNative <- wc(ndat = subset(data_for_fstat, population != "Breeding"), diploid = TRUE)
wc_fst_WildBreeding <- wc(ndat = subset(data_for_fstat, population != "NativeSelection"), diploid = TRUE)
wc_fst_NativeBreeding <- wc(ndat = subset(data_for_fstat, population != "Wild"), diploid = TRUE)


wc_fst_WildNative$per.loc %>%
  rownames_to_column("marker") %>%
  left_join(., snp_info) %>%
  mutate(chrom = parse_number(chrom), even_chrom = chrom %% 2 == 0) %>%
  ggplot(aes(x = pos, y = FST, color = even_chrom)) +
  geom_hline(yintercept = quantile(wc_fst_WildNative$per.loc$FST, 0.999), lty = 2) +
  geom_point() +
  scale_color_discrete(guide = FALSE) +
  scale_x_continuous() +
  facet_grid(~ chrom, scales = "free_x", switch = "x") +
  theme_classic() +
  theme(panel.spacing.x = unit(0, "line"))

wc_fst_WildBreeding$per.loc %>%
  rownames_to_column("marker") %>%
  left_join(., snp_info) %>%
  mutate(chrom = parse_number(chrom), even_chrom = chrom %% 2 == 0) %>%
  ggplot(aes(x = pos, y = FST, color = even_chrom)) +
  geom_hline(yintercept = quantile(wc_fst_WildBreeding$per.loc$FST, 0.999), lty = 2) +
  geom_point() +
  scale_color_discrete(guide = FALSE) +
  scale_x_continuous() +
  facet_grid(~ chrom, scales = "free_x", switch = "x") +
  theme_classic() +
  theme(panel.spacing.x = unit(0, "line"))

wc_fst_NativeBreeding$per.loc %>%
  rownames_to_column("marker") %>%
  left_join(., snp_info) %>%
  mutate(chrom = parse_number(chrom), even_chrom = chrom %% 2 == 0) %>%
  ggplot(aes(x = pos, y = FST, color = even_chrom)) +
  geom_hline(yintercept = quantile(wc_fst_NativeBreeding$per.loc$FST, 0.999, na.rm = TRUE), lty = 2) +
  geom_point() +
  scale_color_discrete(guide = FALSE) +
  scale_x_continuous() +
  facet_grid(~ chrom, scales = "free_x", switch = "x") +
  theme_classic() +
  theme(panel.spacing.x = unit(0, "line"))



# Conduct mantel test between pairwise genetic distance and geogra --------

## Add pairwise genetic distance data?


# Convert the pairwise genetic distance df to a matrix

pairwise_dist <- pairwise_dist_data %>%
  select(contains("individual"), distance) %>%
  spread(individual1, distance) %>%
  as.data.frame() %>%
  column_to_rownames("individual2") %>%
  as.dist()

# Pairwise distance of individuals
pairwise_geo_dist <- germ_meta1 %>%
  filter(individual %in% row.names(geno_mat_wild)) %>%
  select(individual, latitude = origin_latitude, longitude = origin_longitude) %>%
  crossing(., ., .name_repair = tidyr_legacy) %>%
  mutate(gcd = geosphere::distGeo(p1 = select(., longitude, latitude), p2 = select(., longitude1, latitude1))) %>%
  select(individual, individual1, gcd) %>%
  # Filter
  filter_at(vars(contains("individual")), all_vars(. %in% unique(pairwise_dist_data$individual1))) %>%
  spread(individual, gcd) %>%
  as.data.frame() %>%
  column_to_rownames("individual1") %>%
  as.dist()

# Conduct mantel test
mantel_out <- vegan::mantel(xdis = pairwise_dist, ydis = pairwise_geo_dist, )

mantel_out

# It is significant, but the problem is the relationships is not really linear




# Save everything ---------------------------------------------------------



# Save the results
save("global_marker_fst", "global_marker_fst2", "global_fst", "global_fst2", "pairwise_wc_fst", "spa_out",
     "ld_wild_decay", "ld_wild_df", "pairwise_dist_data", "wc_fst_WildNative", "wc_fst_WildBreeding", "wc_fst_NativeBreeding",
     file = file.path(result_dir, "population_genetics_stats.RData"))

