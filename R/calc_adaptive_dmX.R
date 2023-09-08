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
