#' calc_multiadaptive_cov_plot - Calcualte covariance matrix for mult-adaptive covariance plot
#'
#' @param a Rate parameter for OU model
#' @param sigma2_y Variance of Y
#' @param beta Slope parameter
#' @param x X axis value
#' @param Z_adaptive Number of adaptive predictors
#' @param n_reg Number of regimes
#'
#' @return Variance/covariance matrix
#' @export
#'
calc_multiadaptive_cov_plot<-function(a,sigma2_y,beta,x,Z_adaptive,n_reg){
  #Assuming n_reg>1
  ti<-1
  var_opt = sum(beta^2)
  term0 = ((var_opt + sigma2_y) / (2 * a)) * (1 - exp( -2 * a * (1-x/2))) * exp(-a * x)
  term1 = (1 - exp(-a * ti)) / (a * ti)
  term2 = exp(-a * x) * (1 - exp(-a * ti)) / (a * ti)
  Vt<-term0 + var_opt * ((1-x/2) * term1 * term1 - ((1 - exp(-a * (1-x/2))) / a) * (term2 + term2))
  return(Vt)
}
