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
#' @export
part_cor <- function(data, feat, exposure, covars){
  if(length(covars) > 0){
    tmp_data <- data[,c(feat, exposure, covars)]
    tmp_data_dum <- fastDummies::dummy_cols(tmp_data, remove_selected_columns=TRUE)
  } else{
    tmp_data_dum <- data[,c(feat, exposure, covars)]
  }
  tmp_pcor <- ppcor::pcor(tmp_data_dum)
  return(c(feat=feat, coef=tmp_pcor$estimate[1,2], pval=as.numeric(tmp_pcor$p.value[1,2]), type="part_cor"))
}

#' Compute Pearson correlations
#'
#' This function calculates the Pearson correlation between a feature and an exposure, adjusting for covariates.
#'
#' @param data A data frame containing the variables of interest.
#' @param feat A string specifying the feature of interest.
#' @param exposure A string specifying the exposure variable.
#' @param covars A character vector of covariate names to adjust for.
#'
#' @return A named vector with:
#' \itemize{
#'   \item \code{feat}: the feature name,
#'   \item \code{coef}: the estimated coefficient for the exposure variable,
#'   \item \code{pval}: the p-value associated with the exposure variable,
#'   \item \code{type}: the string "pearson".
#' }
#'
#' @details This function fits a linear model of the form: \code{feat ~ exposure + covariates}.
#' The estimated coefficient and p-value for the exposure variable are returned as a measure of
#' correlation, adjusted for covariates.
#'
#' @examples
#' pearson_cor(data = my_data, feat = "metabolite1", exposure = "alcohol", covars = c("age", "sex"))
#'
#' @export
pearson_cor <- function(data, feat, exposure, covars){
  tmp_data <- data[,c(feat, exposure, covars)]
  fmla <- as.formula(paste0(feat, " ~ ", paste(c(exposure, covars),collapse =" + ")))
  res_lm <- lm(formula = fmla, data = tmp_data)
  summary_lm <- summary(res_lm)
  coef <- summary_lm$coefficients[exposure,"Estimate"]
  pval <- summary_lm$coefficients[exposure,"Pr(>|t|)"]
  return(c(feat=feat, coef=as.numeric(coef), pval=as.numeric(pval), type="pearson"))
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
#' @export
filter_cor <- function(cor_res, thresh=0.05, fdr=T){
  if(fdr){ 
    sig_idx <- which(cor_res$fdr <= thresh)
    sig_res <- cor_res$feat[sig_idx]
  } else{
    sig_idx <- which(cor_res$pval <= thresh)
    sig_res <- cor_res$feat[sig_idx]
  }
  return(sig_res)
}

#' Compute and Filter Correlation Features
#'
#' Computes correlations (Pearson or partial) between a set of features and an exposure variable,
#' optionally in parallel. Returns either the full correlation results or a filtered list of significant features.
#'
#' @param data A data frame containing the features, exposure, and covariates.
#' @param features A character vector or numeric vector indicating which features to test. If numeric, column indices are used.
#' @param exposure A character string indicating the exposure variable.
#' @param covars A character vector specifying the covariates to adjust for (only used for partial correlation).
#' @param parallel Logical; whether to run the computation in parallel.
#' @param ncores Integer; number of cores to use for parallel processing (only used if \code{parallel = TRUE}).
#' @param cortype Character; type of correlation to compute. Options are \code{"pearson"} or \code{"partial"}.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Pearson correlation on selected features
#' get_cor_feats(data, features = 1:100, exposure = "alcohol", covars = NULL,
#'               parallel = TRUE, ncores = 4, cortype = "pearson")
#'
#' # Partial correlation with covariates
#' get_cor_feats(data, features = c("feat1", "feat2"), exposure = "alcohol",
#'               covars = c("age", "sex"), parallel = FALSE, cortype = "partial")
#' }
get_cor_feats <- function(data, features, exposure, covars, parallel, ncores=NA, cortype="pearson"){
  if(parallel){
    `%faire%` <- foreach::`%dopar%`
    doParallel::registerDoParallel(ncores)
  } else{`%faire%` <- foreach::`%do%`}
  if(cortype == "pearson"){corfun <- pearson_cor}
  if(cortype == "partial"){corfun <- part_cor}
  if(is.numeric(features)){features <- colnames(data)[features]}
  # Get correlation for each feature
  cor_res <- foreach::foreach(feat=features, .combine='rbind', .packages = c('dplyr', 'ppcor', 'tidyr')) %faire%
    corfun(data, feat, exposure, covars)
  cor_res <- as.data.frame(cor_res)
  cor_res$fdr <- stats::p.adjust(cor_res$pval)
  return(cor_res)
}
