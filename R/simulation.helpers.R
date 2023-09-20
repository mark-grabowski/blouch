#' calc_adaptive_dmX - calculate adaptive predictor matrix for Blouch
#'
#' @param phy Phylogeny in NEXUS format
#' @param a Rate parameter for OU model
#' @param X Predictor
#'
#' @return Adaptive predictor matrix
#' @export
calc_adaptive_dmX<-function(phy,a,X){
  ts<-ts_fxn(phy)
  T_term<-ts[[3]]
  N<-length(T_term);
  if(is.null(dim(X))==FALSE){Z<-dim(X)[2]}else{Z<-1}
  rho<-(1 - (1 - exp(-a * T_term))/(a * T_term))
  rhos<-matrix(rep(rho,Z),nrow=N,ncol=Z)
  dmX<-X * rhos
  return(dmX)
}

#' calc_mixed_dmX - calculate mixed direct effect and adaptive predictor matrix for Blouch
#'
#' @param phy Phylogeny in NEXUS format
#' @param a Rate parameter for OU model
#' @param X Predictor
#' @param Z_direct Number of direct effect predictor traits
#' @param Z_adaptive Number of adaptive predictor traits
#'
#' @return Predictor matrix
#' @export
#'
calc_mixed_dmX<-function(phy,a,X,Z_direct,Z_adaptive){
  ts<-ts_fxn(phy)
  T_term<-ts[[3]]
  N<-length(T_term);
  rho<-(1 - (1 - exp(-a * T_term))/(a * T_term))
  rhos<-matrix(rep(rho,Z_adaptive),nrow=N,ncol=Z_adaptive)
  dmX<-cbind(X[,1:Z_direct],X[,(Z_direct+1):(Z_adaptive+Z_direct)]*rhos)
  return(dmX)
}

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
  N<-dim(ta)[1];
  ti<-matrix(T_term,length(T_term),N);
  var_opt<-sum(as.matrix(beta^2)%*%sigma2_x)
  term0<-(var_opt + sigma2_y) / (2 * a) * (1 - exp( -2 * a * ta)) * exp(-a * tij)
  term1<-(1 - exp(-a * ti)) / (a * ti)
  term2<-exp(-a * tja) * (1 - exp(-a * ti)) / (a * ti)
  Vt<-term0 + var_opt * (ta * term1 * t(term1) - ((1 - exp(-a * ta)) / a) * (term2 + t(term2)))
  return(Vt)
}

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

