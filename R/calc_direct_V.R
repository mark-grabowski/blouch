#' calc_direct_V Calculate direct effect model V/CV matrix
#'
#' @param phy Phylogeny in NEXUS format
#' @param sigma2_y Variance of Y
#' @param a Rate parameter for OU model
#'
#' @return Variance/Covariance marix
#' @export
#'
calc_direct_V<-function(phy, sigma2_y, a){
  ts<-ts_fxn(phy)
  ta<-ts[[1]]
  tij<-ts[[2]] #Same as Cophenetic Distance matrix
  Vt<-sigma2_y /( 2 * a) * ((1 - exp(-2 * a * ta)) * exp(-a * tij)) #ta - time from root to tips, tij  - total time separating spcies
  return(Vt)
}
