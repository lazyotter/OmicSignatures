#' Predict Exposure Signatures from Cross-Validated Model
#'
#' This function uses a cross-validated model object (`cv.glmnet`) to predict exposure signatures for a test set. It returns the predicted values for both the `lambda.min` and `lambda.1se` values from the cross-validation.
#'
#' @param data A data frame or matrix containing the full set of features and exposure variable.
#' @param crossval A `cv.glmnet` object resulting from cross-validation that contains the fitted model.
#' @param test_idx A vector of indices specifying the set for which predictions are to be made.
#' @param features A character vector or column indices specifying the feature names used for prediction.
#' @param exposure A character string specifying the name of the exposure (response variable) in the `data` used for comparison with predictions.
#'
#' @return A list containing two matrices:
#' \itemize{
#'   \item \code{l1se}: A matrix with the actual exposure data and predicted signature using \code{lambda.1se}.
#'   \item \code{lmin}: A matrix with the actual exposure data and predicted signature using \code{lambda.min}.
#' }
#' Each matrix has the exposure values in the first column and the predicted signature values in the second column.
#'
#' @details The function predicts exposure signatures using the fitted `cv.glmnet` model. Predictions are made for both the `lambda.1se` and `lambda.min` values. The returned matrices combine actual exposure values and predicted signatures for comparison.
#' 
#' @examples
#' \dontrun{
#' # Example usage
#' data <- your_data_frame
#' crossval <- cv.glmnet(x = data[train_idx, features], y = data[train_idx, exposure])
#' test_idx <- sample(1:nrow(data), 50)  # Specify test set indices
#' pred_results <- pred_signature(data, crossval, test_idx, features, exposure)
#' }
#'
#' @seealso \code{\link[glmnet]{cv.glmnet}}
#' @export
#'
pred_signature <- function(data, crossval, test_idx, features, exposure) {
  # Ensure that crossval is a cv.glmnet object
  if (!inherits(crossval, "cv.glmnet")) {
    stop("The 'crossval' parameter must be a cv.glmnet object.")
  }
  
  # Extract coefficients
  coefs_l1se <- coef(crossval, s = crossval$lambda.1se)
  coefs_lmin <- coef(crossval, s = crossval$lambda.min)
  
  # Extract feature coefficients only (no covariates)
  features_coefs_l1se <- coefs_l1se[rownames(coefs_l1se) %in% features, , drop = FALSE]
  features_coefs_lmin <- coefs_lmin[rownames(coefs_lmin) %in% features, , drop = FALSE]
  
  # Filter for test set + metabolomic features
  data_tmp <- as.matrix(data[test_idx, rownames(features_coefs_l1se)])
  
  # Get prediction based on lasso coefficients 
  pred_l1se <- data_tmp[,rownames(features_coefs_l1se)] %*% features_coefs_l1se
  pred_lmin <- data_tmp[,rownames(features_coefs_lmin)] %*% features_coefs_lmin
   
  # Get actual data for exposure from test set
  actual_data <- data[test_idx, exposure]
  
  # Combine actual data and predicted signature values
  pred_data <- list(
    l1se = cbind(as.data.frame(actual_data), predsig = as.data.frame(pred_l1se)),
    lmin = cbind(as.data.frame(actual_data), predsig = as.data.frame(pred_lmin))
  )
  
  return(pred_data)
}
