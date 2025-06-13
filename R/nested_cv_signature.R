#' Perform Nested Cross-Validation to Derive and Evaluate a Metabolic Signature
#'
#' This function performs nested cross-validation (CV) to derive a penalized regression signature
#' of an exposure using selected features, and evaluates its performance in held-out test folds.
#' The inner CV is used to tune the model (via `create_signature`), and predictions are evaluated
#' in the outer test fold.
#'
#' @param data A data.frame or matrix with samples in rows and variables (including features and exposure) in columns.
#' @param features A character vector of column names corresponding to the features used to derive the signature.
#' @param exposure A character scalar naming the column of the exposure variable to predict.
#' @param n_folds Integer. Number of folds to use in the outer CV loop (default is 5).
#' @param ... Additional arguments passed to `create_signature`.
#'
#' @return A data.frame containing model performance metrics (Pearson correlation) for each lambda selection method (`l1se` and `lmin`).
#' It includes the following columns:
#' \describe{
#'   \item{pval}{P-value of the Pearson correlation test.}
#'   \item{est}{Estimated Pearson correlation coefficient between predicted signature and true exposure.}
#'   \item{lower_95}{Lower bound of the 95% confidence interval for the correlation.}
#'   \item{upper_95}{Upper bound of the 95% confidence interval for the correlation.}
#'   \item{method}{Lambda selection method used ("l1se" or "lmin").}
#'   \item{exposure}{The name of the exposure variable.}
#'   \item{null_folds}{Number of outer CV folds where no signature was obtained.}
#' }
#'
#' @details
#' This function assumes that signatures are constructed using `create_signature()` and evaluated using `pred_signature()`.
#' If no features survive univariate filtering in a given fold, no signature is computed and a warning is issued.
#'
#' @examples
#' \dontrun{
#' results <- nested_cv_signature(data = my_data, features = metabolite_names, exposure = "alcohol")
#' }
#'
#' @export
nested_cv_signature <- function(data, features, exposure, n_folds=5, ...){
  # Split dataset into n_folds for nested CV
  folds <- split_k_folds(data, n_folds)
  
  pred_cv <- vector("list", 2)
  names(pred_cv) <- c("l1se", "lmin") 
  model_score <- matrix(data = NA, nrow = 0, ncol = 5)
  
  null_folds <- 0
  # Evaluate model for each fold
  for(fold in 1:length(folds)){
    message("(Nested CV) Outer fold ", fold, " out of ", n_folds, " in progress...\n")
    
    # Split data into train & test
    test <- folds[[fold]]
    train <- setdiff(c(1:nrow(data)), test)
    
    # Do 5-fold CV
    train_cv <- create_signature(data = data,
                                 train_idx = train, 
                                 features = features, 
                                 exposure = exposure, 
                                 ...)
    
    if(!is.null(train_cv$crossval)){
      # Get prediction of signatures for each method
      tmp_pred_cv <- pred_signature(data = data,
                                   crossval = train_cv$crossval, 
                                   test_idx = test,
                                   features = features,
                                   exposure = exposure)
      
      # Save results if features selected by lmin or l1se were not none
      for(method in names(tmp_pred_cv)){
        pred_cv[[method]] <- rbind(pred_cv[[method]], tmp_pred_cv[[method]])
        }
    } else{
      null_folds = null_folds + 1
      warning("...no signature obtained for fold ", fold, " out of ", n_folds, ". Interpret performance with caution.")
      }
  }
  
  if(null_folds != n_folds){
    for(method in names(pred_cv)){ # for both l1se and lmin
      cor <- cor_sum(cor.test(x = pred_cv[[method]][,2], # Predicted signature
                              y = pred_cv[[method]][,1])) # Orig value of exposure
      # Add lambda selection method and variable name
      cor[["method"]] <- method
      # Remove auto name of cor
      rownames(cor) <- c()
      
      # Update table of model performance
      model_score <- rbind(model_score, cor)
    }
    model_score$exposure <- exposure
    model_score$null_folds <- null_folds
    message("Nested CV completed. Pearson correlation between exposure & signature in test sets are: \n...for lambda.min: ", round(subset(model_score, method == "lmin")$est,2), "\n...for lambda.1se: ", round(subset(model_score, method == "l1se")$est,2), "\n")
    message("Signature could not be computed in ", null_folds, " out of ", n_folds, " folds.")
    return(model_score)
  } else{
    message("Nested CV was unable to be completed on all folds. Either the signal is weak or there are <= 2 features univariately associated with your exposure.")
    model_score <- data.frame(
      pval = NA,
      est = NA,
      lower_95 = NA,
      upper_95 = NA,
      method = NA,
      exposure = exposure,
      null_folds = null_folds,
      stringsAsFactors = FALSE
    )
    return(model_score)
  }
}