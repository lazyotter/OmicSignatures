% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bootstrap_utils.R
\name{get_sel_freq}
\alias{get_sel_freq}
\title{Get Selection Frequency from Bootstrap Results}
\usage{
get_sel_freq(boot_res, sel_method = "lmin")
}
\arguments{
\item{boot_res}{A list containing bootstrap results, where selection counts are stored in elements named with the prefix \code{"bootstrap_coefs_"} followed by the selection method.}

\item{sel_method}{A string specifying the selection method to use. Options are "lmin" (default) or "l1se".}
}
\value{
A data frame with two columns: \code{count}, the number of times each feature was selected, and \code{sel_freq}, the selection frequency (count divided by the number of bootstrap iterations).
}
\description{
This function calculates the selection frequency of features from bootstrap results.
}
\examples{
boot_res <- list(bootstrap_coefs_lmin = matrix(sample(0:1, 100, replace = TRUE), nrow = 10))
get_sel_freq(boot_res)
}
