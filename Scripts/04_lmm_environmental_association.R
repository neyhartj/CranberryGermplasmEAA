# CranberryGermplasmEAA
#
# Conduct a EAA using the environmental GWAS approach
#

# Load packages
library(sommer)
library(tidyverse)
library(readxl)
library(rrBLUP)
library(snps)
library(qvalue)
library(neyhart)

# What FDR threshold to use?
fdr_thresh <- 0.20


# Read in the base data
load(file.path(data_dir, "population_metadata_and_genotypes.RData"))
# Load the worldclim data
load(file.path(data_dir, "germplasm_origin_bioclim_data.RData"))

# Adjust elevation (zero if missing)
eaa_environmental_data1 <- eaa_environmental_data %>%
  select(location_abbr, all_of(eaa_environmental_vars$variable)) %>%
  mutate(elevation = ifelse(is.na(elevation), 0, elevation))


# Rename genotype matrices
geno_mat_wild <- geno_mat_wild_filter3
geno_hmp_wild <- geno_hmp_wild_filter3


# Prepare bioclim and marker data -----------------------------------------


# Prepare the bioclim data as if it were trait data
germplasm_bioclim_data <- pop_metadata %>%
  filter(category == "Wild", !is.na(latitude)) %>%
  select(individual, latitude, longitude) %>%
  left_join(., eaa_environmental_data1) %>%
  as.data.frame() %>%
  filter(individual %in% colnames(K_wild))


## Prepare the marker data
##
## Subset the individual with bioclim data from the marker matrix
geno_mat_wild1 <- geno_mat_wild[germplasm_bioclim_data$individual,]

# Split the marker matrix by chromosome
geno_mat_list <- geno_hmp_wild %>%
  select(marker, chrom) %>%
  split(.$chrom) %>% map("marker") %>%
  map(~geno_mat_wild1[,.x, drop = FALSE])

# Convert the geno_mat_wild1 matrix into a data.frame for use in the GWAS function
geno <- t(geno_mat_wild1) %>%
  as.data.frame() %>%
  rownames_to_column("marker") %>%
  inner_join(select(snp_info, marker, chrom, pos), .) %>%
  as.data.frame()

# K_wild was calculated using many markers, including those with lower MAF; use this
# relationship matrix for the GWAS
K_wild_use <- K_wild[row.names(geno_mat_wild1), row.names(geno_mat_wild1)]

# Eigenvalues from K
K_eigens <- prcomp(K_wild_use) %>%
  broom::tidy() %>%
  filter(PC %in% 1:3) %>%
  mutate(PC = paste0("PC", PC)) %>%
  spread(PC, value) %>%
  rename(individual = row) %>%
  as.data.frame()



# Run eGWAS ---------------------------------------------------------------



# For each variable, fit the following models:
# 1. K
# 2. K + PCs
# 3. K + latitude
# 4. K + PCs + latitude
# 5. K + latitude + longitude
# 6. K + PCs + latitude + longitude

# Create a tibble of model parameters to check
egwas_models_df <- tribble(
  ~model, ~K, ~PCs, ~covariates,
  "model1", TRUE, 0, character(),
  "model2", TRUE, 3, character(),
  "model3", TRUE, 0, c("latitude"),
  "model4", TRUE, 3, c("latitude"),
  "model5", TRUE, 0, c("latitude", "longitude"),
  "model6", TRUE, 3, c("latitude", "longitude"),
) %>%
  mutate(npar = map2_dbl(PCs, covariates, ~1 + 1 + .x + length(.y) + 1),
         output = list(NULL))


# Subset to create the pheno data
pheno <- germplasm_bioclim_data %>%
  select(individual, all_of(eaa_environmental_vars$variable)) %>%
  as.data.frame()


# Iterate over the models in the list
for (i in seq_len(nrow(egwas_models_df))) {

  # Pull out other terms
  covariates <- egwas_models_df$covariates[[i]]
  nPC <- egwas_models_df$PCs[i]

  # Empty list for scores
  scores_list <- list()

  # Iterate over the variables
  for (biovar in eaa_environmental_vars$variable) {

    # If biovar is in "covariates," remove it
    covariates1 <- setdiff(covariates, biovar)

    # Subset pheno
    pheno1 <- pheno[,c("individual", biovar, covariates1)]
    # Set covariates1 to NULL if empty
    if (length(covariates1) == 0) covariates1 <- NULL

    # Run the GWAS
    gwas_out <- GWAS(pheno = pheno1, geno = geno, n.PC = nPC, min.MAF = 0, P3D = TRUE,
                     plot = FALSE, fixed = covariates1, K = K_wild_use)

    # Add the scores to a matrix
    scores_list[[biovar]] <- gwas_out

  }

  # Merge the list
  scores_df <- reduce(scores_list, full_join, by = c("marker", "chrom", "pos"))

  # Add this data to the output df
  egwas_models_df$output[[i]] <- scores_df

} # Close the loop



## Combine the scores for different models
gwas_out <- egwas_models_df %>%
  select(-K, -PCs, -covariates) %>%
  unnest(output) %>%
  gather(variable, score, eaa_environmental_vars$variable) %>%
  # Calculate p-values and uniform p-values for comparison
  mutate(p_value = 10^-score) %>%
  split(list(.$variable, .$model)) %>%
  map_dfr(~{
    .x2 <- filter(.x, score != 0) %>%
      mutate(unif_score = qqscore(score)) %>%
      select(marker, chrom, pos, model, variable, score, p_value, unif_score)

    left_join(.x, .x2, by = c("model", "marker", "chrom", "pos", "variable", "score", "p_value"))
  })

# P-value inflation factor
p_lambda <- gwas_out %>%
  filter(!is.na(p_value)) %>%
  group_by(variable, model) %>%
  summarize(p_value_inflation = QCEWAS::P_lambda(p = p_value), .groups = "drop") %>%
  # Assign the model with inflation factor closes to 1 as the best model
  group_by(variable) %>%
  mutate(best_model = model[which.min(abs(1 - p_value_inflation))],
         best_model = ifelse(model == best_model, "*", "")) %>%
  ungroup()


## Examine Q-Q plots

# Plot
qq_plots <- gwas_out %>%
  left_join(., p_lambda) %>%
  # filter(variable %in% eaa_environmental_vars$variable[1:4]) %>%
  mutate(best_model = best_model == "*") %>%
  ggplot(aes(x = unif_score, y = score, color = model, shape = best_model)) +
  geom_abline(slope = 1, intercept = 0, lwd = 0.5) +
  geom_point(size = 0.75) +
  facet_wrap(~ variable) +
  theme_bw()

# Save
ggsave(filename = "egwas_qqplots.jpg", plot = qq_plots, path = fig_dir,
       height = 25, width = 10, dpi = 100)


# Select significant SNPs
eaa_gwas_sigmar <- gwas_out %>%
  filter(model == "model2") %>%
  split(.$variable) %>%
  map_df(~filter(.x, p_value <= sommer:::fdr(p_value, fdr.level = fdr_thresh)$fdr.10))



# For each variable, use the population-structure corrected climate information to
# fit a model that estimates SNP R2
mlm_out_list <- list()

for (vari in unique(eaa_gwas_sigmar$variable)) {

  # Get the number of PCs and covariates
  terms <- subset(egwas_models_df, model == "model2")
  nPCs <- terms$PCs
  covariates <- terms$covariates[[1]]

  # Create a new pheno object
  pheno1 <- pheno %>%
    left_join(select(K_eigens, individual, PC1, PC2, PC3), by = c("individual", "PC1", "PC2", "PC3"))


  # Create the model formula
  fixed <- reformulate(c(paste0(rep("PC", nPCs), seq_len(nPCs)), covariates), response = vari)
  # Build matrices to fit the mixed model
  mf <- model.frame(formula = modelr::add_predictors(fixed, ~ individual), data = pheno1)

  y <- model.response(mf)
  X <- model.matrix(fixed, mf)
  Z <- model.matrix(~ -1 + individual, mf)
  K1 <- K_wild_use[mf$individual, mf$individual]

  # Fit the mixed model
  mm_out <- mixed.solve(y = y, Z = Z, K = K1, X = X, method = "REML", SE = TRUE)

  # Calculate adjusted variable values
  adj_y <- y - ((X %*% mm_out$beta) + (Z %*% mm_out$u))

  # Regress these values on the significant SNPs

  # Subset markers with p-values below this score
  sig_mar <- subset(eaa_gwas_sigmar, variable == vari, marker, drop = TRUE)

  # Create a model frame with the SNP genotypes
  mf2 <- as.data.frame(cbind(y = y, y_adj = as.numeric(adj_y), geno_mat_wild1[,sig_mar, drop = FALSE]))
  # Run the regression
  mr_out <- lm(reformulate(sig_mar, response = "y_adj", intercept = FALSE), mf2)
  # Stepwise removal of markers
  mr_out_step <- step(mr_out, trace = 0)
  # Do the same thing, but regress the original y values
  mr_out_yorig <- lm(reformulate(sig_mar, response = "y", intercept = FALSE), mf2)
  mr_out_step_yorig <- step(mr_out_yorig, trace = 0)

  # Calculate r squared for retained markers
  mr_out_step_adjRsquared <- summary(mr_out_step)$adj.r.square
  mr_out_step_yorig_adjRsquared <- summary(mr_out_step_yorig)$adj.r.square

  # Calculate r squared for PCs and polygenic background
  pop_str_r2 <- as.numeric(var((X %*% mm_out$beta))) / var(y)
  polygen_r2 <- (mm_out$Vu / var(y))
  resid_r2 <- 1 - (pop_str_r2 + polygen_r2)

  # Calculate the variance explained by each component of the model
  mlm_variance_explained <- tibble(
    term = c("population_structure", "polygenic_background", "significant_markers"),
    nLevels = as.numeric(c(NA, NA, length(coef(mr_out_step)))),
    total_variance_explained = c(pop_str_r2, polygen_r2, mr_out_step_adjRsquared * resid_r2),
    resid_variance_explained = as.numeric(c(NA, NA, mr_out_step_adjRsquared)),
    nLevels_noLMM = as.numeric(c(NA, NA, length(coef(mr_out_step_yorig)))),
    total_variance_explained_noLMM = as.numeric(c(NA, NA, mr_out_step_yorig_adjRsquared))
  )

  # Add everything to the list
  mlm_out_list[[vari]] <- mlm_variance_explained

}

# Combine
mlm_out_df <- imap_dfr(mlm_out_list, ~mutate(.x, variable = .y)) %>%
  select(variable, names(.))


# ## Save the results
# save("gwas_out", "p_lambda", "mlmm_out_list", file = file.path(result_dir, "eaa_gwas_results.RData"))

# Edit egwas_models_df
egwas_models_df <- egwas_models_df %>%
  select(-output)


save("gwas_out", "p_lambda", "eaa_gwas_sigmar", "egwas_models_df", "mlm_out_df",
     file = file.path(result_dir, "eaa_gwas_results.RData"))










