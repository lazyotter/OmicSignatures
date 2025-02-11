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

#' Compute Partial Correlation
#'
#' This function calculates the partial correlation between a feature and an exposure, adjusting for covariates.
#'
#' @param data A data frame containing the variables of interest.
#' @param feat A string specifying the feature of interest.
#' @param exposure A string specifying the exposure variable.
#' @param covars A character vector of covariate names to adjust for.
#'
#' @return A named vector with the feature name, partial correlation coefficient, p-value, and correlation type.
#'
#' @importFrom dplyr select
#' @importFrom fastDummies dummy_cols
#' @importFrom ppcor pcor
#' @export
part_cor <- function(data, feat, exposure, covars){
  tmp_data <- crc_RC_all %>% dplyr::select(all_of(c(feat, exposure, covars)))
  tmp_data_dum <- dummy_cols(tmp_data, remove_selected_columns=TRUE)
  tmp_pcor <- pcor(tmp_data_dum)
  return(c(feat=feat, coef=tmp_pcor$estimate[1,2], pval=tmp_pcor$p.value[1,2], type="part_cor"))
}

#' Filter Significant Correlations
#'
#' This function filters correlation results based on a significance threshold, optionally applying FDR correction.
#'
#' @param cor_res A data frame containing correlation results with p-values.
#' @param thresh A numeric threshold for significance (default: 0.05).
#' @param fdr A logical indicating whether to apply FDR correction (default: TRUE).
#'
#' @return A character vector of significant feature names.
#'
#' @importFrom dplyr mutate filter
#' @importFrom stats p.adjust
#' @export
filter_cor <- function(cor_res, thresh=0.05, fdr=T){
  if(fdr){ 
    sig_feats <- cor_res %>% as.data.frame %>% mutate(fdr=p.adjust(pval)) %>%
      filter(fdr <= thresh) %>% pull(feat) 
  } else{
    sig_feats <- as.data.frame(cor_res) %>% filter(pval <= thresh)
  }
}

#' Compute Partial Correlation for Multiple Features
#'
#' This function computes partial correlations for multiple features and filters significant ones.
#'
#' @param data A data frame containing the variables of interest.
#' @param features A character vector of feature names.
#' @param exposure A string specifying the exposure variable.
#' @param covars A character vector of covariate names to adjust for.
#' @param parallel A logical indicating whether to run computations in parallel.
#' @param ncores An integer specifying the number of cores to use (default: NA).
#'
#' @return A filtered data frame of significant features after partial correlation analysis.
#'
#' @importFrom foreach %dopar% %do% foreach
#' @importFrom doParallel registerDoParallel
#' @export
get_pcor_feats <- function(data, features, exposure, covars, parallel, ncores=NA){
  if(parallel){
    `%faire%` <- `%dopar%`
    registerDoParallel(ncores)
  } else{`%faire%` <- `%do%`}
  if(is.numeric(features)){features <- colnames(data)[features]}
  cor_res <- foreach(feat=features, .combine='rbind', .packages = c('dplyr', 'ppcor', 'tidyr')) %faire%
    part_cor(data, feat, exposure, covars)
  feats_filt <- filter_cor(cor_res)
  return(feats_filt)
}
