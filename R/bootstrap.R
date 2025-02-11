#' Perform Bootstrapping for Signature Creation
#'
#' This function performs bootstrapping to evaluate the stability of feature selection and signature performance 
#' in a cross-validation setting using LASSO or a similar approach.
#'
#' @param data A data frame or matrix containing the dataset. Rows represent samples, and columns represent features, covariates, and the exposure variable.
#' @param features A character vector of column names in `data` representing the features to include in the model.
#' @param exposure A character string representing the column name in `data` corresponding to the exposure variable.
#' @param n_boot An integer specifying the number of bootstrap iterations.
#' @param covars A character vector of column names in `data` representing additional covariates to include in the model. Defaults to an empty vector.
#' @param folds An integer specifying the number of folds for cross-validation. Defaults to 5.
#' @param parallel A logical value indicating whether to use parallel processing. Defaults to `FALSE`.
#' @param cores An integer specifying the number of CPU cores to use if `parallel = TRUE`. Defaults to 18.
#'
#' @return A list with the following components:
#' \describe{
#'   \item{bootstrap_coefs_l1se}{A sparse matrix containing bootstrap counts for each feature selected using the 1-SE rule.}
#'   \item{bootstrap_coefs_lmin}{A sparse matrix containing bootstrap counts for each feature selected using the minimum lambda.}
#'   \item{bootstrap_cors}{A data frame containing the number of selected features and the correlations between predicted signatures and the exposure variable for each bootstrap iteration.}
#' }
#'
#' @examples
#' \dontrun{
#' # Example dataset
#' data <- data.frame(matrix(rnorm(1000), nrow = 100, ncol = 10))
#' colnames(data) <- c(paste0("feature", 1:8), "exposure", "covar")
#' results <- bootstrap_signature(
#'   data = data,
#'   features = paste0("feature", 1:8),
#'   exposure = "exposure",
#'   n_boot = 100,
#'   covars = "covar",
#'   folds = 5,
#'   parallel = TRUE,
#'   cores = 4
#' )
#' }
#'
#' @export
bootstrap_signature <- function(data, features, exposure, n_boot, covars = c(), folds = 5, parallel = F, cores = 18){
  
  # Create a dataframe with counts of each feature for each iteration
  boot_feat_count.l1se <- Matrix::Matrix(0, nrow = n_boot, ncol = length(features), dimnames = list(c(1:n_boot), features), sparse = T)
  boot_feat_count.lmin <- Matrix::Matrix(0, nrow = n_boot, ncol = length(features), dimnames = list(c(1:n_boot), features), sparse = T)
  
  # Create an empty data frame with column names
  results_df <- data.frame(
    num_feat_l1se = numeric(),
    num_feat_lmin = numeric(),
    cor_l1se = numeric(),
    pval_l1se = numeric(),
    cor_lmin = numeric(),
    pval_lmin = numeric(),
    stringsAsFactors = FALSE
  )
  
  t1 <- Sys.time()
  for(i in 1:n_boot){
    set.seed(i)
    # Sample with replacement same number of samples as original data
    index_bootstrap <- sample(1:nrow(data), nrow(data), replace = T)
    # Get count of each sample
    index_freq_table <- table(index_bootstrap) %>% as.data.frame()
    
    # Index of rows selected from bootstrap
    rows_bootstrap <- index_freq_table$index_bootstrap %>% as.character() %>% as.numeric()
    
    # Create matrix of features, covars, and exposure
    tmp_data <- data[rows_bootstrap, c(features, covars, exposure)]
    # Run cross-validation
    if(parallel){registerDoMC(cores = cores)}

    tmp_res <- create_signature(data = tmp_data, 
                                train_idx = c(1:nrow(tmp_data)), 
                                features = features, 
                                exposure = exposure,
                                covars = covars, 
                                folds = folds,
                                parallel = parallel, 
                                seed = i,
                                cores = cores,
                                weights = as.matrix(index_freq_table$Freq))
    # tmp_res <- cv.glmnet(x = data_x, y = data_y, weights = index_freq_table$Freq, lambda.min.ratio = 0.001,
    #                      type.measure = "mse",
    #                      nlambda = 100, parallel = parallel, nfolds = folds)
    tmp_coef_l1se <- coef(tmp_res$crossval, s="lambda.1se")
    tmp_coef_lmin <- coef(tmp_res$crossval, s="lambda.min")
    
    cols <- which(rownames(tmp_coef_l1se) %in% features) # Saving coefs for features only
    boot_feat_count.l1se[i,] <- tmp_coef_l1se[cols] # Exclude intercept
    boot_feat_count.lmin[i,] <- tmp_coef_lmin[cols] # Exclude intercept

    # Get prediction of signatures on all data
    tmp_pred_cv <- pred_signature(data = tmp_data,
                                 crossval = tmp_res$crossval,
                                 test_idx = 1:nrow(tmp_data),
                                 features = features,
                                 exposure = exposure)

    # Get correlation of signature to original exposure
    if(sum(is.na(tmp_pred_cv$l1se)) > 0){
      tmp_cor_l1se <- data.frame(estimate = NA, p.value = NA)
    } else {
      tmp_cor_l1se <- cor.test(tmp_pred_cv$l1se[,1], tmp_pred_cv$l1se[,2])}
    if(sum(is.na(tmp_pred_cv$lmin)) > 0){
      tmp_cor_lmin <- data.frame(estimate = NA, p.value = NA)
    } else {
      tmp_cor_lmin <- cor.test(tmp_pred_cv$lmin[,1], tmp_pred_cv$lmin[,2])}
    
    # Save correlation results, number of features selected in this iteration
    ## Number of features selected with l1se or lmin
    print(which(tmp_coef_l1se[cols] != 0))
    num_feat_l1se <- length(which(tmp_coef_l1se[cols] != 0)) - 1 # Exclude intercept
    num_feat_lmin <- length(which(tmp_coef_lmin[cols] != 0)) - 1 # Exclude intercept
    ## Append values to the data frame
    results_df <- rbind(results_df, data.frame(
      num_feat_l1se = num_feat_l1se,
      num_feat_lmin = num_feat_lmin,
      cor_l1se = tmp_cor_l1se$estimate,
      pval_l1se = tmp_cor_l1se$p.value,
      cor_lmin = tmp_cor_lmin$estimate,
      pval_lmin = tmp_cor_lmin$p.value
    ))
  }
  t2 <- Sys.time()
  cat("Bootstrapping completed in", difftime(t2,t1), "\n")

  res <- list(bootstrap_coefs_l1se = boot_feat_count.l1se,
              bootstrap_coefs_lmin = boot_feat_count.lmin,
              bootstrap_cors = results_df)
  return(res)
}
