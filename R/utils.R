#' Summarize Correlation Test Results
#'
#' This function extracts key statistics from a correlation test object and returns them in a tidy data frame format.
#'
#' @param cor A correlation test object, typically the output of \code{cor.test()} in base R. 
#' This object should contain p-values, correlation estimates, and confidence intervals.
#'
#' @return A data frame with the following columns:
#' \itemize{
#'   \item \code{pval}: The p-value of the correlation test.
#'   \item \code{est}: The correlation estimate.
#'   \item \code{lower_95}: The lower bound of the 95% confidence interval for the correlation estimate.
#'   \item \code{upper_95}: The upper bound of the 95% confidence interval for the correlation estimate.
#' }
#'
#' @details 
#' This function is a simple utility for summarizing the results of \code{cor.test()} in a structured format 
#' that is easy to use for downstream analysis or reporting.
#'
#' @examples
#' \dontrun{
#' # Perform a correlation test
#' cor_test_result <- cor.test(mtcars$mpg, mtcars$wt)
#' 
#' # Summarize the correlation test results
#' summary <- cor_sum(cor_test_result)
#' print(summary)
#' }
#' 
#' @export
#'
cor_sum <- function(cor){
  pval <- cor$p.value
  est <- cor$estimate
  lower_95 <- cor$conf.int[1]
  upper_95 <- cor$conf.int[2]
  return(data.frame(pval = pval, est = est, lower_95 = lower_95, upper_95 = upper_95))
}

#' Split Data into K Folds
#'
#' This function splits a dataset into \code{n_folds} approximately equal-sized folds for cross-validation. 
#' Each fold contains a unique set of indices, ensuring no overlap.
#'
#' @param data A data frame or matrix to be split into folds. The function uses the number of rows to create the folds.
#' @param n_folds Integer specifying the number of folds to create.
#'
#' @return A list of length \code{n_folds}, where each element is a vector of row indices corresponding to one fold.
#'
#' @details 
#' This function randomly splits the data into \code{n_folds}. The size of each fold is approximately \code{floor(nrow(data) / n_folds)}. 
#' The last fold contains any remaining rows due to rounding.
#'
#' @examples
#' \dontrun{
#' # Split a dataset into 5 folds
#' folds <- split_k_folds(iris, 5)
#' print(folds)
#' }
#' 
#' @export
#'
split_k_folds <- function(data, n_folds){
  # sample(nrow(data))
  res <- list()
  leftover_index <- c(1:nrow(data))
  size_fold <- floor(nrow(data)/n_folds)
  
  # Get indices for each fold of size_fold
  for(i in 1:(n_folds-1)){
    fold <- sample(x = leftover_index, size = size_fold, replace = F) %>% sort()
    res[[i]] <- fold
    leftover_index <- setdiff(leftover_index, fold)
  }
  # Set last fold to all indices left
  res[[i+1]] <- leftover_index
  
  return(res)
}

#' Convert Covariates to Factors
#'
#' This function converts specified columns in a data frame to factors.
#'
#' @param data A data frame containing the columns to be converted.
#' @param covars A character vector or numeric vector specifying the columns to convert to factors.
#'
#' @return A data frame with the specified columns converted to factors.
#' @examples
#' # Example with character vector of column names
#' data <- data.frame(a = c(1, 2, 3), b = c("x", "y", "z"), c = c(TRUE, FALSE, TRUE))
#' covar_to_factor(data, c("a", "c"))
#'
#' # Example with numeric indices
#' data <- data.frame(a = c(1, 2, 3), b = c("x", "y", "z"), c = c(TRUE, FALSE, TRUE))
#' covar_to_factor(data, c(1, 3))
#'
#' @export
#'
covar_to_factor <- function(data, covars){
  data[, covars] <- lapply(data[, covars], as.factor)
  return(data)
}

#' Convert Covariates to Numeric
#'
#' This function converts specified columns in a data frame to numeric.
#'
#' @param data A data frame containing the columns to be converted.
#' @param covars A character vector or numeric vector specifying the columns to convert to numeric.
#'
#' @return A data frame with the specified columns converted to numeric.
#'
#' @examples
#' # Example with character vector of column names
#' data <- data.frame(a = c("1", "2", "3"), b = c("4", "5", "6"), c = c(TRUE, FALSE, TRUE))
#' covar_to_numeric(data, c("a", "b"))
#'
#' # Example with numeric indices
#' data <- data.frame(a = c("1", "2", "3"), b = c("4", "5", "6"), c = c(TRUE, FALSE, TRUE))
#' covar_to_numeric(data, c(1, 2))
#'
#' @export
#'
covar_to_numeric <- function(data, covars){
  data[, covars] <- lapply(data[, covars], as.numeric)
  return(data)
}