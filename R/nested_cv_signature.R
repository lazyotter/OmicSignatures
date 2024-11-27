#' Nested Cross-Validation for Signature Evaluation
#'
#' This function performs nested cross-validation to evaluate the performance of a LASSO-derived omics signature. 
#' It uses a specified exposure and feature set, returning a summary of model performance.
#'
#' @param data A data frame containing the features (predictor variables) and exposure (response variable).
#' @param features A character vector of column names specifying the features to be included in the model.
#' @param exposure A character string specifying the name of the exposure variable (response) in `data`.
#'
#' @return A data frame summarizing model performance, with columns for the method (`lambda.min` or `lambda.1se`), 
#' estimated correlation, p-value, and 95% confidence interval.
#'
#' @details 
#' The function performs the following steps:
#' \enumerate{
#'   \item Splits the dataset into 5 folds for nested cross-validation using \code{\link{split_k_folds}}.
#'   \item For each fold:
#'     \itemize{
#'       \item Splits the data into training and testing subsets.
#'       \item Fits a LASSO model on the training set using \code{\link{create_signature}}.
#'       \item Predicts the exposure using the fitted model on the test set via \code{\link{pred_signature}}.
#'       \item Saves predictions for both \code{lambda.min} and \code{lambda.1se}.
#'     }
#'   \item Calculates correlation metrics between predicted and actual exposure values for each method 
#'         using \code{\link{cor_sum}} and \code{\link[stats]{cor.test}}.
#'   \item Returns a performance summary with method name, correlation estimate, p-value, and confidence intervals.
#' }
#'
#' @examples
#' \dontrun{
#' # Example dataset
#' data <- my_data_frame
#' features <- colnames(data)[2:100] # Example feature names
#' exposure <- "my_exposure"
#' 
#' # Perform nested cross-validation
#' results <- nested_cv_signature(data, features, exposure)
#' print(results)
#' }
#'
#' @seealso 
#' \code{\link{split_k_folds}}, 
#' \code{\link{create_signature}}, 
#' \code{\link{pred_signature}}, 
#' \code{\link{cor_sum}}
#'
#' @export
nested_cv_signature <- function(data, features, exposure, n_folds=5, ...){
  # Split dataset into n_folds for nested CV
  folds <- split_k_folds(data, n_folds)
  
  pred_cv <- vector("list", 2)
  names(pred_cv) <- c("l1se", "lmin") 
  model_score <- matrix(data = NA, nrow = 0, ncol = 5)
  
  # Evaluate model for each fold
  for(fold in 1:length(folds)){
    cat("Outer fold", fold, "out of", n_folds, "in progress...\n")
    
    # Split data into train & test
    test <- folds[[fold]]
    train <- setdiff(c(1:nrow(data)), test)
    
    # Do 5-fold CV
    train_cv <- create_signature(data = data,
                                 train_idx = train, 
                                 features = features, 
                                 exposure = exposure, 
                                 ...)
    
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
    }
  
  
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
return(model_score)
}