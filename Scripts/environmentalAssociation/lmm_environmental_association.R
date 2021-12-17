# Germplasm collection environmental association
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

# Load the startup script
source("startup.R")


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




# # Use LRTs to determine which model to run
#
# full_models <- pheno %>%
#   gather(variable, value, -individual) %>%
#   group_by(variable) %>%
#   do({
#     df <- .
#     df1 <- left_join(df, K_eigens, by =  "individual") %>%
#       left_join(., pheno[,c("individual", "latitude", "longitude")], by = "individual")
#
#     # Copy egwas_models_df
#     egwas_models_df1 <- egwas_models_df %>%
#       select(-output) %>%
#       mutate(logLik = as.numeric(NA))
#
#     # Iterate over models
#     for (i in seq_len(nrow(egwas_models_df1))) {
#       nPC <- egwas_models_df1$PCs[i]
#       covariates <- setdiff(egwas_models_df1$covariates[[i]], unique(df1$variable))
#       random_terms <- "individual"
#       fixed_terms <- c("1", covariates, paste0(rep("PC", nPC), seq_len(nPC)))
#
#       # Create a model frame
#       mf <- model.frame(reformulate(c(random_terms, fixed_terms), response = "value"),
#                         data = df1)
#       # Model matrices
#       Z <- model.matrix(~ -1 + individual, mf)
#       X <- model.matrix(reformulate(fixed_terms), mf)
#
#       # Fit the model
#       fit <- mixed.solve(y = model.response(mf), Z = Z, K = K_wild_use, X = X, method = "ML")
#       # Return the LL and df
#       egwas_models_df1$logLik[i] <- fit$LL
#
#     }
#
#     # Calculate p-values
#     mod2vmod1 <- pchisq(q = 2 * (egwas_models_df1$logLik[2] - egwas_models_df1$logLik[1]),
#                         df = egwas_models_df1$npar[2] - egwas_models_df1$npar[1], lower.tail = FALSE)
#     mod3vmod1 <- pchisq(q = 2 * (egwas_models_df1$logLik[3] - egwas_models_df1$logLik[1]),
#                         df = egwas_models_df1$npar[3] - egwas_models_df1$npar[1], lower.tail = FALSE)
#     mod4vmod2 <- pchisq(q = 2 * (egwas_models_df1$logLik[4] - egwas_models_df1$logLik[2]),
#                         df = egwas_models_df1$npar[4] - egwas_models_df1$npar[2], lower.tail = FALSE)
#     mod5vmod3 <- pchisq(q = 2 * (egwas_models_df1$logLik[5] - egwas_models_df1$logLik[3]),
#                         df = egwas_models_df1$npar[5] - egwas_models_df1$npar[3], lower.tail = FALSE)
#     mod6vmod4 <- pchisq(q = 2 * (egwas_models_df1$logLik[6] - egwas_models_df1$logLik[4]),
#                         df = egwas_models_df1$npar[6] - egwas_models_df1$npar[4], lower.tail = FALSE)
#
#     egwas_models_df1 %>%
#       mutate(p_value = c(as.numeric(NA), mod2vmod1, mod3vmod1, mod4vmod2, mod5vmod3, mod6vmod4))
#
#   }) %>% ungroup()





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


#
# # For each variable, fit a multi-locus model to calculate SNP effects
# mlmm_out_list <- list()
#
#
#
# for (vari in unique(gwas_out$variable)) {
#   # Subset the scores for the best model
#   best_model <- subset(p_lambda, variable == vari & best_model == "*", model, drop = TRUE)
#
#   scores1 <- gwas_out %>%
#     filter(variable == vari, model == best_model)
#
#   # Get the number of PCs and covariates
#   terms <- subset(egwas_models_df, model == best_model)
#   nPCs <- terms$PCs
#   covariates <- terms$covariates[[1]]
#
#   # Calculate the FDR p-value
#   fdr_p <- sommer:::fdr(p = scores1$p_value, fdr.level = fdr_thresh)$fdr.10
#   fdr_p <- ifelse(fdr_p > 1, 0, fdr_p)
#
#   # Subset markers with p-values below this score
#   sig_mar <- subset(scores1, p_value <= fdr_p, marker, drop = TRUE)
#
#   if (length(sig_mar) == 0) next
#
#   # Create a new pheno object
#   pheno1 <- pheno %>%
#     left_join(select(pc_eigenvec, individual, PC1, PC2, PC3)) %>%
#     left_join(., rownames_to_column(as.data.frame(geno_mat_wild1[,sig_mar, drop = FALSE]), "individual"))
#
#   # Create the model formula
#   fixed <- as.formula(paste(vari, "~ 1")) %>%
#     modelr::add_predictors(f = ., as.formula(paste0("~", paste0(c(paste0(rep("PC", nPCs), seq_len(nPCs)), covariates, sig_mar), collapse = "+"))))
#
#   # Build matrices to fit the mixed model
#   mf <- model.frame(formula = modelr::add_predictors(fixed, ~ individual), data = pheno1)
#
#   y <- model.response(mf)
#   X <- model.matrix(fixed, mf)
#   Z <- model.matrix(~ -1 + individual, mf)
#   K1 <- K_wild_use[mf$individual, mf$individual]
#
#   # Fit the mixed model
#   mlmm_out <- mixed.solve(y = y, Z = Z, K = K1, X = X, method = "REML", SE = TRUE)
#
#   u <- mlmm_out$u
#   which_x <- which(! colnames(X) %in% colnames(geno_mat_wild1))
#   b <- X[,which_x,drop = FALSE] %*% mlmm_out$beta[which_x]
#
#   # Calculate marginal y (y - polygenic effect)
#   y_marginal <- y - as.vector(b) - u
#
#   # Test SNP effects
#   fixef <- cbind(beta = mlmm_out$beta, se = mlmm_out$beta.SE)
#   # p-values
#   p_vals <- pchisq(q = fixef[,"beta"]^2 / fixef[,"se"]^2, df = 1, lower.tail = FALSE)
#   fixef <- cbind(fixef, p_value = p_vals)
#
#   ## While loop
#   non_sig <- any(fixef[sig_mar, "p_value"] > mlmm_p_thresh)
#   while(non_sig) {
#     # Exclude the least significant marker; rerun the model
#     sig_mar <- setdiff(sig_mar, names(which.max(fixef[row.names(fixef) %in% scores1$marker,"p_value"])[1]))
#     X <- X[,setdiff(colnames(X), intersect(setdiff(colnames(X), sig_mar), colnames(geno_mat_wild1))), drop = FALSE]
#
#     # Fit the mixed model
#     mlmm_out <- mixed.solve(y = y, Z = Z, K = K1, X = X, method = "REML", SE = TRUE)
#
#     # Test SNP effects
#     fixef <- cbind(beta = mlmm_out$beta, se = mlmm_out$beta.SE)
#     # p-values
#     p_vals <- pchisq(q = fixef[,"beta"]^2 / fixef[,"se"]^2, df = 1, lower.tail = FALSE)
#     fixef <- cbind(fixef, p_value = p_vals)
#
#     non_sig <- any(fixef[sig_mar, "p_value"] > mlmm_p_thresh)
#
#
#   }
#
#
#   # Calculate R2 as correlation of predicted with observed
#   beta_mat <- matrix(data = fixef[,"beta"], nrow = nrow(X), ncol = ncol(X), byrow = TRUE,
#                      dimnames = list(NULL, row.names(fixef)))
#   # Individual markers
#   r2_cor <- apply(X = X * beta_mat, MARGIN = 2, FUN = function(y_hat) cor(y_hat, y)^2)
#
#   # All markers
#   y_hat_all_mar <- X[,sig_mar, drop = FALSE] %*% fixef[,"beta"][sig_mar]
#   r2_all_mar <- cor(y_hat_all_mar, y)^2
#
#   # Also calculate R2 for the random polygenic background
#   r2_polygenic <- cor(mlmm_out$u, y)^2
#
#
#   # Effects
#   effect_plot_df <- rownames_to_column(as.data.frame(cbind(fixef, r2_cor)), "term")
#   r2_additional <- tibble(term = c("significant_markers", "all_markers"), r2_cor = c(r2_all_mar, r2_polygenic))
#
#   # Add everything to the list
#   mlmm_out_list[[vari]] <- list(effect_plot_df, r2_additional,
#                                 y_marginal = rownames_to_column(as.data.frame(y_marginal), "individual"))
#
#
# }
#
#
# # Count associations
# mlmm_out_df <- mlmm_out_list %>%
#   map(1) %>% # Get the first item
#   map(~filter(.x, str_detect(term, "^S"))) %>%
#   imap_dfr(~mutate(.x, variable = .y))
#
# # Bioclim variables
# mlmm_out_df %>%
#   group_by(variable) %>%
#   count() %>%
#   as.data.frame()
#
# # Markers
# mlmm_out_df %>%
#   group_by(term) %>%
#   count() %>%
#   as.data.frame()
#
# ## Marker related to temperature
# mlmm_out_df %>%
#   left_join(., eaa_environmental_vars) %>%
#   group_by(class, term) %>%
#   count() %>%
#   as.data.frame()



# ## Save the results
# save("gwas_out", "p_lambda", "mlmm_out_list", file = file.path(result_dir, "eaa_gwas_results.RData"))

# Edit egwas_models_df
egwas_models_df <- egwas_models_df %>%
  select(-output)


save("gwas_out", "p_lambda", "egwas_models_df", file = file.path(result_dir, "eaa_gwas_results.RData"))




