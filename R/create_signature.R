#' Format cv.glmnet Object from create_signature
#'
#' Extracts selected features from a `cv.glmnet` object using the lambda values for minimum error (`lambda.min`) and one standard error (`lambda.1se`).
#'
#' @param df_crossval A `cv.glmnet` object, typically obtained from `create_signature`.
#'
#' @return A list with two elements:
#' \itemize{
#'   \item \code{crossval}: The original `cv.glmnet` object.
#'   \item \code{results}: A list with selected variables for `lambda.min` and `lambda.1se`.
#' }
#' 
#'
format_cvglmnet <- function(df_crossval, features){
  
  resultsCV_tmp <- list()
  
  resultsCV_tmp$lmin <- extract_features(df_crossval, "lambda.min", features)
  
  resultsCV_tmp$l1se <- extract_features(df_crossval, "lambda.1se", features)
  
  res <- list(df_crossval, resultsCV_tmp)
  
  names(res) <- c("crossval", "results")
  return(res)
}

#' Extract Selected Features from a cv.glmnet Object at Specified Lambda
#'
#' This function extracts the selected features from a `cv.glmnet` object based on the coefficients at a specified lambda value. 
#' It returns the names of the features that have non-zero coefficients, excluding the intercept.
#'
#' @param df_crossval A `cv.glmnet` object, typically the result of running `cv.glmnet` on a dataset. This object contains the cross-validation results and the fitted model.
#' @param lambda A character string specifying the lambda value to use for feature selection. This can be either "lambda.min" or "lambda.1se" to select features based on the respective lambda value.
#' @param features A character vector of column names specifying all feature names.
#'
#' @return A character vector containing the names of the features selected at the specified lambda, excluding the intercept.
#'
#' @details This function uses the `cv.glmnet` object to extract the coefficients at the specified lambda value. It identifies the features with non-zero coefficients and returns their names. The intercept is excluded from the result.
#'
#' @examples
#' \dontrun{
#' # Example usage of extract_features
#' df_crossval <- glmnet::cv.glmnet(x = training_data, y = response_data)
#' selected_features_min <- extract_features(df_crossval, "lambda.min")
#' selected_features_1se <- extract_features(df_crossval, "lambda.1se")
#' print(selected_features_min)
#' print(selected_features_1se)
#' }
#'
#' @seealso \code{\link[glmnet]{cv.glmnet}}, \code{\link{format_cvglmnet}}
#' @export
#'
extract_features <- function(df_crossval, lambda, features){
  coef_res <- predict(df_crossval$glmnet.fit, type = "coefficients", s = df_crossval[[lambda]])
  
  # Extract non-zero feature names at lambda
  sig_features <- names(which(coef_res[, 1] != 0)[-1])  # Exclude intercept
  # Select features only (if there are covariates)
  sig_features <- sig_features[which(sig_features %in% features)]
  return(sig_features)
}

#' Create Omics-Based Exposure Signature Using Cross-Validated LASSO
#'
#' Fits a cross-validated LASSO model using `cv.glmnet` to identify a signature of features associated with a specified outcome. Optionally supports parallel processing for faster computation.
#'
#' @param data A data frame containing both features and the outcome variable.
#' @param train_idx A vector of row indices specifying the training subset of `data`.
#' @param features A character vector of column names specifying the features (predictor variables) to be included in the model.
#' @param exposure A character string specifying the name of the exposure/outcome variable (response) in `data`.
#' @param covars A character vector specifying covariates to be included in the model. These will not be penalized.'
#' @param filter Logical value indicating whether to filter features used in the model for those partially correlated with the exposure.
#' @param folds Integer specifying the number of folds for cross-validation. Default is 5.
#' @param parallel Logical value indicating whether to use parallel processing. Default is \code{FALSE}.
#' @param seed Integer specifying the seed for set.seed. Default is 1.
#' @param cores Integer specifying the number of cores to use if \code{parallel = TRUE}. Default is 18.
#'
#' @return A list with two elements:
#' \itemize{
#'   \item \code{crossval}: The fitted `cv.glmnet` object.
#'   \item \code{results}: A list of selected variables based on `lambda.min` and `lambda.1se`.
#' }
#' 
#' @details This function performs LASSO regression using cross-validation to select optimal features. If \code{parallel = TRUE}, the function uses \code{doMC} for parallel computation with the specified number of cores.
#'
#' @examples
#' \dontrun{
#' data <- my_data_frame
#' train_index <- sample(1:nrow(data), size = floor(0.8 * nrow(data)))
#' feature_names <- colnames(data)[2:100] # Example feature names
#' outcome_name <- "my_outcome"
#' signature <- create_signature(data, train_index, feature_names, outcome_name, folds = 10, parallel = TRUE, n_cores = 4)
#' }
#'
#' @seealso \code{\link[glmnet]{cv.glmnet}}, \code{\link[doMC]{registerDoMC}}, \code{\link{format_cvglmnet}}
#' @export
#'
create_signature <- function(data, train_idx, features, exposure, covars = c(), filter = F, folds = 5, parallel = F, cores = 18, seed = 1, ...){
  t1 <- Sys.time()
  n_covars <- 0
  # Filtering features for significant partial correlation w/exposure
  if(filter){
    feats_final <- get_pcor_feats(data, features, exposure, covars, parallel, cores)
  } else{
    feats_final <- features
  }
  
  # Feature/covariate dataframe
  if(length(covars) > 0){
    data_covars <- fastDummies::dummy_cols(data[, covars],
                                           remove_selected_columns = T, 
                                           remove_first_dummy = T)
    n_covars <- ncol(data_covars)
    data_tmp <- cbind(data[train_idx, feats_final], data_covars[train_idx,]) %>% as.matrix()
  } else {
    data_tmp <- as.matrix(data[train_idx, feats_final])
  }
  # Create penalty vector
  penalty_factors <- c(rep(1, length(feats_final)), rep(0, n_covars))
  # Exposure dataframe
  outcome_tmp <- data[train_idx, exposure] %>% as.matrix()
  colnames(outcome_tmp) <- exposure
  
  # Start parallel instance if using parallelization
  if(parallel){
    doMC::registerDoMC(cores = cores)
  }
  set.seed(seed)
  crossval_tmp <- glmnet::cv.glmnet(x = data_tmp, y = outcome_tmp, 
                                    lambda.min.ratio = 0.001,
                                    type.measure = "mse", nlambda = 100, 
                                    parallel = parallel, nfolds = folds,
                                    penalty.factor = penalty_factors, ...)
  if(is.numeric(features)){
    features <- colnames(data[,features])
  }
  res <- format_cvglmnet(crossval_tmp, features)
  
  t2 <- Sys.time()
  cat("Cross-validated signature completed in", difftime(t2,t1), "\n")
  
  return(res)
}
