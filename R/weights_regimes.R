#' weights_regimes - for internal Blouch use
#' For individual lineage, sum up the segments in each regimes
#' @param a Rate parameter from OU model
#' @param lineage Individual regime values for lineage
#'
#' @return Return named vector with regimes weights for individual lineage
#' @export
#'
weights_regimes <- function(a, lineage) {#For individual lineage, sum up the segments in each regimes
  #nt <- lineage$nodes_time
  res <- weights_segments(a, lineage) ## Rcpp wrapper, equivalent to above commented out code
  w <- vapply(lineage$which.regimes, function(e) sum(e*res), FUN.VALUE = 0) ## Sum coefficients for which the respective regime is equal
  return(w) #Return named vector with regimes weights for individual lineage
}
