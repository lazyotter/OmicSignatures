#' Get Selection Frequency from Bootstrap Results
#'
#' This function calculates the selection frequency of features from bootstrap results.
#'
#' @param boot_res A list containing bootstrap results, where selection counts are stored in elements named with the prefix `"bootstrap_coefs_"` followed by the selection method.
#' @param sel_method A string specifying the selection method to use. Options are "lmin" (default) or "l1se".
#'
#' @return A data frame with two columns: `count`, the number of times each feature was selected, and `sel_freq`, the selection frequency (count divided by the number of bootstrap iterations).
#'
#' @export
#'
#' @examples
#' boot_res <- list(bootstrap_coefs_lmin = matrix(sample(0:1, 100, replace = TRUE), nrow = 10))
#' get_sel_freq(boot_res)
get_sel_freq <- function(boot_res, sel_method="lmin"){
  sel_name <- paste0("bootstrap_coefs_", sel_method)
  count <- colSums(boot_res[[sel_name]] > 0)
  count <- as.data.frame(count) %>% dplyr::arrange(desc(count))
  count <- count %>% dplyr::mutate(sel_freq = count/nrow(boot_res[[sel_name]]))
  return(count)
}