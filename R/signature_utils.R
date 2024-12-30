#' Extract Non-Zero Coefficients from a LASSO Model
#'
#' This function extracts the non-zero coefficients from a LASSO model fit stored 
#' within a signature object. It returns the coefficients along with the feature names 
#' in a data frame format, excluding the intercept.
#'
#' @param alc_signature A list object containing the cross-validation model fit, 
#' typically created using the `create_signature` function. The cross-validation model 
#' should be stored in the `crossval` element.
#' @param lambda Character or numeric. The penalty parameter (lambda) at which to 
#' extract the coefficients. Default is `"lambda.min"` for the lambda that gives 
#' the minimum cross-validation error.
#'
#' @return A data frame with two columns:  
#' \describe{
#'   \item{`Feature`}{The name of the feature (omics variable).}
#'   \item{`Coefficient`}{The estimated coefficient for the feature.}
#' }
#' The intercept is excluded from the returned data frame.
#'
#' @examples
#' # Example usage with a fitted signature object:
#' # Assuming `alc_signature` is an object created using `create_signature()`
#' signature <- get_lasso_coefficients(alc_signature, lambda = "lambda.min")
#' print(signature)
#'
#' @export
get_lasso_coefficients <- function(alc_signature, lambda = "lambda.min") {
  tmp_model <- alc_signature$crossval
  # Get coefficients for the specified lambda
  coefficients <- coef(tmp_model, s = lambda)
  
  # Convert to a data frame
  coef_df <- as.data.frame(as.matrix(coefficients))
  coef_df <- coef_df[coef_df[, 1] != 0, , drop = FALSE]
  colnames(coef_df) <- "Coefficient"
  coef_df$Feature <- rownames(coef_df)
  rownames(coef_df) <- NULL
  coef_df <- coef_df[-1,] # Remove intercept
  return(coef_df)
}
