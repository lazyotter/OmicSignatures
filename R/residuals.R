#' Calculate Residuals for Features Adjusted for Covariates and Exposure
#'
#' This function computes residuals for a given set of features after adjusting 
#' for specified covariates and an exposure variable. It uses linear models for adjustment 
#' and supports parallel processing to improve performance when analyzing multiple features.
#'
#' @param data A data frame containing the features, covariates, and exposure variables.
#' @param features A character vector of feature names (e.g., metabolite names) or a numeric vector of column indices for which residuals are to be calculated.
#' @param covars A character vector of covariate column names or a numeric vector of column indices to adjust for in the linear models.
#' @param exposure A character string specifying the exposure variable column name or a numeric index of the exposure variable to include in the models.
#' @param id_cols A character vector of column names or a numeric vector of column indices to use as identifiers (e.g., sample IDs) in the resulting data frame.
#' @param ctrl_idx A numeric vector specifying the row indices of control samples. If empty, all rows are treated as controls.
#' @param parallel A logical value indicating whether to use parallel processing. Defaults to `FALSE`.
#' @param cores An integer specifying the number of cores to use for parallel processing. Defaults to 18.
#' 
#' @return A data frame containing the residuals for each feature, along with the identifier columns and remaining data columns.
#'
#' @details
#' The function fits a linear model for each feature using the formula:
#' \code{`feature_name ~ exposure + covariates`}.
#' Residuals are computed for control samples (if specified) and optionally for case samples.
#' Residuals are adjusted by adding back the fixed effect of the exposure variable.
#' 
#' If `ctrl_idx` is provided, the linear model is fit only on the control samples, and residuals computed for the case samples using the control model.
#' 
#'
#' @examples
#' # Example dataset
#' data <- data.frame(
#'   ID = 1:100,
#'   feature1 = rnorm(100),
#'   feature2 = rnorm(100),
#'   covar1 = runif(100),
#'   covar2 = runif(100),
#'   exposure = sample(c(0, 1), 100, replace = TRUE)
#' )
#' 
#' # Specify parameters
#' features <- c("feature1", "feature2")
#' covars <- c("covar1", "covar2")
#' exposure <- "exposure"
#' id_cols <- "ID"
#' 
#' # Compute residuals
#' residuals <- get_residuals(
#'   data = data,
#'   features = features,
#'   covars = covars,
#'   exposure = exposure,
#'   id_cols = id_cols,
#'   parallel = FALSE
#' )
#'
#' @import dplyr
#' @import foreach
#' @import doParallel
#' @import lme4
#' @export
get_residuals <- function(data, features, covars, exposure, id_cols, ctrl_idx=c(), parallel=F, cores=18){
  
  if (any(is.na(data[, covars]))) {
    warning("There are NA values in the covariates.")
  } 
  
  # Create "blank formula": ~ fixed + random effects
  fmla_lm_base <- paste0("~ ", paste(c(exposure, covars), collapse = " + "))
  
  if(length(ctrl_idx) != 0){
    data_model <- data[ctrl_idx,]
  } else{
    data_model <- data
  }
  
  cat(paste0("Calculating residuals for ", exposure, "...\n"))
  # For each feature, residuals with lm (not using mixed models because not correcting for plate or batch)
  if(parallel){
    # Register parallel backend if needed
    cl <- parallel::makeCluster(cores)
    doParallel::registerDoParallel(cl)
    
  } 
  
  `%>%` <- dplyr::`%>%`

  resid_tmp = foreach(ii=1:length(features), .packages = c("lme4", "tidyverse"), .export = "features", .combine='cbind') %dopar% {
    cat("Computing residuals on", colnames(data[,features])[ii], "...\n")
    
    namefeat_tmp <- features[ii]

    # Formula for each feature (eg. `180.0625@2.8240602` ~ "MMed_Score" + "Alc_Re")
    fmla_lm <- paste(paste0("`",namefeat_tmp,"`"), fmla_lm_base)
    
    # Run linear model on controls
    tmp_res_ctrl <- lm(as.formula(fmla_lm), data=data_model)
    resid_ctrl  <- residuals(tmp_res_ctrl)
    exposure_coef <- coefficients(tmp_res_ctrl)[exposure]
    # Add back the fixed effect into the residual
    res_ctrl <- exposure_coef * data_model %>% dplyr::select(dplyr::all_of(exposure)) + resid_ctrl
    
    # If there is a case/control split, compute residuals on cases
    if(length(ctrl_idx) > 0){
      data_case <- data[-ctrl_idx,]
      tmp_res_case <- predict(tmp_res_ctrl, newdata = data_case)
      resid_case <- data_case[[namefeat_tmp]] - tmp_res_case
      res_case <- exposure_coef * data_case %>% dplyr::select(dplyr::all_of(exposure)) + resid_case
      res <- rbind(res_ctrl, res_case)
    } else{
      res <- res_ctrl
    }
   
    return(res)
  }
  
  if(parallel){parallel::stopCluster(cl)}

  if(length(ctrl_idx) > 0){
    ctrl_names <- data[ctrl_idx,] %>% dplyr::select(dplyr::all_of(id_cols))
    case_names <- data[-ctrl_idx,] %>% dplyr::select(dplyr::all_of(id_cols))
    names <- rbind(ctrl_names, case_names)
    resid <- dplyr::bind_cols(names, resid_tmp)
    # Merge with original covariates
    final <- merge(data %>% dplyr::select(!dplyr::all_of(features)), by = id_cols, no.dups = TRUE, resid)
  } else {
    colnames(resid_tmp) <- features
    final <- dplyr::bind_cols(data %>% dplyr::select(!dplyr::all_of(features)), resid_tmp)
  }
  

  return(final)

  
}