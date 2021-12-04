# Germplasm collection environmental association
#
# Analyze environmental association results
#

# Load packages
library(neyhart)
library(tidyverse)
library(readxl)
library(snps)
library(rrBLUP)

# Load the startup script
source("startup.R")


# Read in the marker data
load(file.path(data_dir, "gc_marker_data.RData"))
load(file.path(data_dir, "germplasm_metadata.RData"))

# Load the worldclim data
load(file.path(data_dir, "germplasm_origin_bioclim_data.RData"))

# Center and scale
eaa_environmental_data1 <- eaa_environmental_data %>%
  select(location_abbr, elevation, latitude, longitude, all_of(eaa_environmental_vars$variable)) %>%
  # mutate_at(vars(contains("bio")), scale, scale = FALSE) %>%
  # mutate_at(vars(contains("bio")), as.numeric) %>%
  mutate(elevation = ifelse(is.na(elevation), 0, elevation))

# Read in variety origin information
native_sel_meta <- read_csv(file = file.path(data_dir, "native_selection_metadata.csv")) %>%
  rename(state = original_planting_state)


# Load the eGWAS results
load(file.path(result_dir, "eaa_gwas_results.RData"))


# Prepare the bioclim data ------------------------------------------------

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


# PCA of wild individuals for determining PCs
K_pca_wild <- prcomp(x = K_wild)

# Gather the eigenvector
pc_eigenvec <- K_pca_wild$rotation %>%
  as.data.frame() %>%
  rownames_to_column("individual") %>%
  as_tibble()

germplasm_bioclim_data1 <- germplasm_bioclim_data %>%
  left_join(., select(pc_eigenvec, individual, PC1, PC2, PC3)) %>%
  left_join(., rownames_to_column(as.data.frame(geno_mat_wild), "individual")) %>%
  rowid_to_column() %>%
  as_tibble()



# Predict environmental conditions based on marker profile ----------------


## First use markers from the mlmm


# Create resamples of the "phenotype" data based on population
lopo_resamples <- germplasm_bioclim_data1 %>%
  split(.$location_abbr) %>%
  map(~list(test = .x, train = filter(germplasm_bioclim_data1, !rowid %in% .x$rowid)))

# what is the size of the smallest training set
min_train_size <- min(sapply(map(lopo_resamples, "train"), nrow))


# Define a function for prediction
gs_predict <- function(train, test, fixed, rand) {

  # Build matrices to fit the mixed model
  mf <- model.frame(formula = modelr::add_predictors(fixed, ~ individual), data = train)
  mf_test <- model.frame(formula = modelr::add_predictors(fixed, ~ individual), data = test)


  # Prepare incidence matrices
  y_train <- model.response(mf)
  X_train <- model.matrix(fixed, mf)
  Z <- model.matrix(rand, mf)
  K1 <- K_wild[levels(mf$individual), levels(mf$individual)]

  # Fit the model
  model_out <- mixed.solve(y = y_train, Z = Z, K = K1, X = X_train, method = "REML", SE = TRUE)

  # Predict
  y_test <- model.response(mf_test)
  X_test <- model.matrix(fixed, mf_test)

  y_test_hat <- (X_test %*% model_out$beta) + as.matrix(model_out$u[mf_test$individual])

  return(y_test_hat)

}






# Create a data.frame to store results
prediction_list1 <- list()

# Iterate over the variables
for (i in seq_along(mlmm_out_list)) {

  # Grab the parameters for the model for this variable
  vari <- names(mlmm_out_list)[i]
  param <- mlmm_out_list[[i]][[1]] %>%
    filter(term != "(Intercept)")

  # Create the model formula
  fixed <- as.formula(paste(vari, "~", paste0(c(1, param$term), collapse = " + ")))
  rand <- ~ -1 + individual

  # Iterate over the training/testing populations
  train_test_out <- list()

  for (j in seq_along(lopo_resamples)) {

    # Separate training and testing sets
    test <- lopo_resamples[[j]]$test
    train <- lopo_resamples[[j]]$train
    train <- mutate(train, individual = factor(individual, levels = c(train$individual, test$individual)))

    y_test_hat <- gs_predict(train = train, test = test, fixed = fixed, rand = rand)

    # Return the testing set and its predictions to the list
    train_test_out[[j]] <- cbind(select(test, rowid, individual, location_abbr, observed_bioclim = all_of(vari)), prediction = y_test_hat[,1])

  } # Close the loop

  # Row bind the train/test
  prediction_results <- bind_rows(map(train_test_out, as_tibble)) %>%
    mutate(squared_error = (observed_bioclim - prediction)^2)

  # Summarize
  prediction_results_summary <- prediction_results %>%
    summarize(cor = cor(observed_bioclim, prediction), rmse = sqrt( mean( squared_error ) ))

  # Add this to the prediction list
  prediction_list1[[vari]] <- list(predictions = prediction_results, summary = prediction_results_summary)

}


## Combine and plot

# Combine summaries
prediction_list2 <- prediction_list1 %>%
  transpose() %>%
  map(., ~imap_dfr(.x = .x, ~mutate(.x, variable = .y)))


g_prediction1 <- prediction_list2$predictions %>%
  ggplot(aes(x = prediction, y = observed_bioclim)) +
  geom_point(aes(color = location_abbr), size = 1) +
  geom_abline(slope = 1, intercept = 0, lwd = 0.5) +
  geom_text(data = mutate(prediction_list2$summary, annotation = paste0("r==", format_numbers(cor, 3))),
            aes(label = annotation, x = Inf, y = -Inf), size = 2, hjust = 1, vjust = -1, parse = TRUE) +
  scale_color_discrete(guide = FALSE) +
  facet_wrap(~ variable, scales = "free") +
  theme_genetics()

# Save
ggsave(filename = "logo_predictions_mlmm_parameters.jpg", plot = g_prediction1, path = fig_dir,
       height = 8, width = 10, dpi = 1000)




## Next look at different p-value thresholds

# Vector of p-value thresholds
p_thresholds <- c(0, 0.001, 0.01, 0.05, 0.10, 0.20, 0.50)

# Create a data.frame to store results
pvalue_thresh_predictions <- crossing(variable = unique(gwas_out$variable), pvalue_threshold = p_thresholds) %>%
  mutate(out = list(NULL))

# Iterate over rows
pb <- progress::progress_bar$new(total = nrow(pvalue_thresh_predictions))
for (i in seq_len(nrow(pvalue_thresh_predictions))) {

  row <- pvalue_thresh_predictions[i,]

  vari <- row$variable
  # Subset the gwas results for markers that meet the threshold
  gwas_out_i <- filter(gwas_out, variable == vari) %>%
    arrange(p_value)
  fdr_thresh <- sommer:::fdr(p = gwas_out_i$p_value, fdr.level = row$pvalue_threshold)$fdr.10
  # Subset the parameters
  param <- tibble(term = head(subset(gwas_out_i, p_value <= fdr_thresh, marker, drop = TRUE), min_train_size - 1))

  # Create the model formula
  fixed <- as.formula(paste(vari, "~", paste0(c(1, param$term), collapse = " + ")))
  rand <- ~ -1 + individual

  # Iterate over the training/testing populations
  train_test_out <- list()

  for (j in seq_along(lopo_resamples)) {

    # Separate training and testing sets
    test <- lopo_resamples[[j]]$test
    train <- lopo_resamples[[j]]$train
    train <- mutate(train, individual = factor(individual, levels = c(train$individual, test$individual)))

    y_test_hat <- gs_predict(train = train, test = test, fixed = fixed, rand = rand)

    # Return the testing set and its predictions to the list
    train_test_out[[j]] <- cbind(select(test, rowid, individual, location_abbr, observed_bioclim = all_of(vari)), prediction = y_test_hat[,1])

  } # Close the loop

  # Row bind the train/test
  prediction_results <- bind_rows(map(train_test_out, as_tibble)) %>%
    mutate(squared_error = (observed_bioclim - prediction)^2)

  # Summarize
  prediction_results_summary <- prediction_results %>%
    summarize(cor = cor(observed_bioclim, prediction), rmse = sqrt( mean( squared_error ) ))

  # Return the results
  pvalue_thresh_predictions$out[[i]] <- mutate(prediction_results_summary, predictions = list(prediction_results))

  pb$tick()
}

pvalue_thresh_predictions1 <- unnest(pvalue_thresh_predictions) %>%
  mutate(pvalue_threshold_fct = fct_inseq(as.character(pvalue_threshold)))

# Plot accuracy by p-value threshold
g_prediction2 <- pvalue_thresh_predictions1 %>%
  ggplot(aes(x = pvalue_threshold, y = cor)) +
  geom_point() +
  geom_line() +
  scale_y_continuous(name = "Prediction accuracy", breaks = pretty) +
  facet_wrap(~variable) +
  labs(subtitle = "Prediction accuracy with different p-value thresholds") +
  theme_genetics()

# Save this
ggsave(filename = "logo_predictions_pvalue_threshold_tests.jpg", plot = g_prediction2, path = fig_dir,
       height = 8, width = 10, dpi = 1000)


# For each variable, print two graphs: accuracy using all markers and accuracy
# using all markers + 0.5 pvalue threshold
pvalue_thresh_predictions2 <- pvalue_thresh_predictions1 %>%
  unnest(predictions)




# Calculate genomewide marker effects and examine germplasm archit --------



# Calculate marker effects
bioclim_marker_eff <- eaa_environmental_vars %>%
  select(variable) %>%
  filter(variable %in% names(germplasm_bioclim_data1)) %>%
  group_by(variable) %>%
  do(marker_effect = {
    row <- .
    vari <- row$variable[[1]]

    # Model formula
    form <- reformulate(c("individual", str_subset(names(germplasm_bioclim_data1), "^PC")), response = vari)
    mf <- model.frame(formula = form, data = germplasm_bioclim_data1)

    # Model matrices
    y <- model.response(mf)
    X1 <- model.matrix(~ 1, mf)
    X2 <- model.matrix(formula(drop.terms(terms(form), 1)), mf)
    Z <- geno_mat_wild[mf$individual,,drop = FALSE]
    # Fit
    fit1 <- mixed.solve(y = y, Z = Z, X = X1)
    fit2 <- mixed.solve(y = y, Z = Z, X = X2)
    lr <- 2 * (fit2$LL - fit1$LL)

    # Return the marker effects and the model
    tibble(marker = names(fit1$u), effect1 = fit1$u, effect2 = fit2$u)

  }) %>% ungroup()


# Calculate chromosome effects per individual per bioclime
bioclim_chromosome_eff <- bioclim_marker_eff %>%
  group_by(variable) %>%
  do(chrom_eff = {
    row <- .
    # Get the marker effect mat
    a <- row$marker_effect[[1]] %>%
      select(marker, effect = effect2) %>%
      as.data.frame() %>%
      column_to_rownames("marker") %>%
      as.matrix()
    a1 <- matrix(t(a), nrow = nrow(geno_mat_wild), ncol = ncol(geno_mat_wild), byrow = TRUE,
                 dimnames = dimnames(geno_mat_wild))
    # Convert to df
    a1_df <- (geno_mat_wild * a1) %>%
      as.data.frame() %>%
      rownames_to_column("individual") %>%
      gather(marker, effect, -individual) %>%
      # Add chromosome information
      left_join(., select(snp_info, marker, chrom), by = "marker") %>%
      group_by(individual, chrom) %>%
      summarize(effect = sum(effect), .groups = "drop") %>%
      # Add population origin
      left_join(., distinct(germplasm_bioclim_data, individual, location_abbr), by = "individual")

    a1_df

  }) %>% ungroup()


# Heatmaps
vari <- "phh2o_topsoil"
df <- subset(bioclim_chromosome_eff, variable == vari, chrom_eff, drop = TRUE)[[1]]

df %>%
  ggplot(aes(x = individual, y = chrom, fill = effect)) +
  geom_tile() +
  scale_fill_gradient2() +
  facet_grid(~ location_abbr, scales = "free_x", space = "free_x") +
  labs(subtitle = subset(eaa_environmental_vars, variable == vari, full_name, drop = TRUE))




