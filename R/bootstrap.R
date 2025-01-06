bootstrap_signature <- function(data, features, exposure, n_boot, covars = c(), folds = 5, parallel = F, cores = 18)
  
  # Create a dataframe with counts of each feature for each iteration
  boot_feat_count.l1se <- Matrix::Matrix(0, nrow = n_boot, ncol = length(features), dimnames = list(features, c(1:n_boot)), sparse = T)
  boot_feat_count.lmin <- Matrix::Matrix(0, nrow = n_boot, ncol = length(features), dimnames = list(features, c(1:n_boot)), sparse = T)
  
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
  for(i in 1:n_bootstrap){
    set.seed(i)
    # Sample with replacement same number of samples as original data
    index_bootstrap <- sample(1:nrow(data), nrow(data), replace = T)
    # Get count of each sample
    index_freq_table <- table(index_bootstrap) %>% as.data.frame()
    
    # Index of rows selected from bootstrap
    rows_bootstrap <- index_freq_table$index_bootstrap %>% as.character() %>% as.numeric()
    
    # Create matrix of features x samples
    data_x <- data[rows_bootstrap, data[,features]] %>% as.matrix()
    # Create matrix of output variable
    data_y <- data[rows_bootstrap, col_name] %>% as.matrix()
    
    # Run cross-validation
    registerDoMC(cores = 32)
    tmp_res <- cv.glmnet(x = data_x, y = data_y, weights = index_freq_table$Freq, lambda.min.ratio = 0.001,
                         type.measure = "mse",
                         nlambda = 100, parallel = T, nfolds = 5)
    tmp_res <- format_cvglmnet(tmp_res)
    
    # Get prediction of signatures on all data
    tmp_pred_cv <- pred_crossval(data = tmp_data,
                                 crossval = tmp_res$crossval,
                                 test_index = 1:nrow(tmp_data),
                                 var = col_name)
    
    # Get correlation of signature to original exposure
    if(sum(is.na(tmp_pred_cv$l1se)) > 0){
      tmp_cor_l1se <- data.frame(estimate = NA, p.value = NA)
    } else {
      tmp_cor_l1se <- cor.test(tmp_pred_cv$l1se[,1], tmp_pred_cv$l1se[,2])}
    if(sum(is.na(tmp_pred_cv$lmin)) > 0){
      tmp_cor_lmin <- data.frame(estimate = NA, p.value = NA)
    } else {
      tmp_cor_lmin <- cor.test(tmp_pred_cv$lmin[,1], tmp_pred_cv$lmin[,2])}
    
    # Save selected features, fit for this iteration
    boot_feat_count.l1se[name, tmp_res$results$l1se] <- boot_feat_count.l1se[name, tmp_res$results$l1se] + 1
    boot_feat_count.lmin[name, tmp_res$results$lmin] <- boot_feat_count.lmin[name, tmp_res$results$lmin] + 1
    
    # Save correlation results, number of features selected in this iteration
    ## Number of features selected with l1se or lmin
    num_feat_l1se <- length(tmp_res$results$l1se)
    num_feat_lmin <- length(tmp_res$results$lmin)
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
  
  assign(paste0(name, "_boot"), results_df)
  t2 <- Sys.time()
  difftime(t1, t2)
  cat("\n")
}
res <- grep(pattern = "boot", objects(), value = "T")