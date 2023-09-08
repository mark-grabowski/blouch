#' calc_adaptive_V - Calculate adaptive variance/covariance matrix
#'
#' @param phy phylogney in NEXUS format
#' @param a Rate parameter for OU model
#' @param sigma2_y Variance of Y
#' @param beta slope
#' @param sigma2_x Brownian-motion parameterof X
#' @param Z_adaptive Number of adaptive predictors
#'
#' @return Variance/covariance matrix
#' @export
#'
calc_adaptive_V<-function(phy, a, sigma2_y, beta,  sigma2_x, Z_adaptive){
  ts<-ts_fxn(phy)
  ta<-ts[[1]]
  tij<-ts[[2]]
  T_term<-ts[[3]]
  tja<-ts[[4]]
  N<-length(T_term);
  ti<-matrix(T_term,length(T_term),N);
  var_opt<-sum(as.matrix(beta^2)%*%sigma2_x)
  term0<-(var_opt + sigma2_y) / (2 * a) * (1 - exp( -2 * a * ta)) * exp(-a * tij)
  term1<-(1 - exp(-a * ti)) / (a * ti)
  term2<-exp(-a * tja) * (1 - exp(-a * ti)) / (a * ti)
  Vt<-term0 + var_opt * (ta * term1 * t(term1) - ((1 - exp(-a * ta)) / a) * (term2 + t(term2)))
  return(Vt)
}
