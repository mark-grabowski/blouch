#' calc_multiadaptive_cov_plot - Calcualte covariance matrix for mult-adaptive covariance plot
#'
#' @param a Rate parameter for OU model
#' @param sigma2_y Variance of Y
#' @param beta Slope parameter
#' @param x X axis value
#'
#' @return Variance/covariance matrix
#' @export
#'
calc_multiadaptive_cov_plot<-function(a,sigma2_y,beta,x){
  ta<-1 #Time from base of phylogeny to MRCA
  ti<-ta+x #Time from base of phylogeny to species i
  tia<-x #Time from species i to the MRCA of two species
  tij<-2*x #Time separating the two species i, j
  var_opt = sum(beta^2)
  term0 = ((var_opt + sigma2_y) / (2 * a)) * (1 - exp( -2 * a * (1-ta/2))) * exp(-a * tij)
  term1 = (1 - exp(-a * ti)) / (a * ti)
  term2 = exp(-a * ti) * (1 - exp(-a * ti)) / (a * ti)
  Vt<-term0 + var_opt * ((1-ta/2) * term1 * term1 - ((1 - exp(-a * (1-ta/2))) / a) * (term2 + term2))
  return(Vt)
}
#' calc_adaptive_cov_plot - Calculate covariance matrix for adaptive covariance plot
#'
#' @param a Rate parameter for OU model
#' @param sigma2_y Variance of Y
#' @param beta Slope parameter
#' @param x X axis value
#'
#' @return Variance/covariance matrix
#' @export
#'
calc_adaptive_cov_plot<-function(a,sigma2_y,beta,x){
  ta<-1 #Time from base of phylogeny to MRCA
  ti<-ta+x #Time from base of phylogeny to species i
  tia<-x #Time from species i to the MRCA of two species
  tij<-2*x #Time separating the two species i, j
  var_opt = sum(beta^2)
  term0 = ((var_opt + sigma2_y) / (2 * a)) * (1 - exp( -2 * a * (1-ta/2))) * exp(-a * tij)
  term1 = (1 - exp(-a * ti)) / (a * ti)
  term2 = exp(-a * ti) * (1 - exp(-a * ti)) / (a * ti)
  Vt<-term0 + var_opt * ((1-ta/2) * term1 * term1 - ((1 - exp(-a * (1-ta/2))) / a) * (term2 + term2))
  return(Vt)
}

