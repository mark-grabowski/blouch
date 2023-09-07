#' calc_adaptive_V - Calculate adaptive variance/covariance matrix
#'
#' @param phy phylogney in phytools format
#' @param a Rate parameter for OU model
#' @param sigma2_y Variance of Y
#' @param ta Tree info
#' @param tij Tree info
#' @param tja Tree info
#' @param T_term Tree info
#' @param beta slope
#' @param sigma2_x Brownian-motion parameterof X
#' @param Z_adaptive Number of adaptive predictors
#' @param n_reg Number of regimes
#'
#' @return Variance/covariance matrix
#' @export
#'
calc_adaptive_V<-function(phy, a, sigma2_y, ta,  tij,  tja,  T_term,  beta,  sigma2_x, Z_adaptive, n_reg){
  N<-dim(ta)[1];
  ts<-ts_fxn(phy)
  ta<-ts[[1]]
  tij<-ts[[2]] #Same as Cophenetic Distance matrix
  T_term<-ts[[3]]
  tja<-ts[[4]]
  ones<-matrix(rep(1,Z_adaptive),nrow=Z_adaptive,ncol=1)
  #beta<-matrix(beta,nrow=n_reg,ncol=1)
  ti<-matrix(T_term,length(T_term),N);
  #if(Z_adaptive>1){
  #  var_opt<-t(beta)%*%sigma2_x%*%ones
  #}
  #else{
  var_opt<-sum(as.matrix(beta^2)%*%sigma2_x)#}
  term0<-(var_opt + sigma2_y) / (2 * a) * (1 - exp( -2 * a * ta)) * exp(-a * tij)
  term1<-(1 - exp(-a * ti)) / (a * ti)
  term2<-exp(-a * tja) * (1 - exp(-a * ti)) / (a * ti)
  Vt<-term0 + var_opt * (ta * term1 * t(term1) - ((1 - exp(-a * ta)) / a) * (term2 + t(term2)))
  return(Vt)
}
