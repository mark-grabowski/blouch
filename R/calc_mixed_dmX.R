#' calc_mixed_dmX - calculate mixed direct effect and adaptive predictor matrix for Blouch
#'
#' @param a Rate parameter for OU model
#' @param T_term T_term parameter for phylogeny
#' @param X Predictor
#' @param Z_direct Number of direct effect predictor traits
#' @param Z_adaptive Number of adaptive predictor traits
#'
#' @return Predictor matrix
#' @export
#'
calc_mixed_dmX<-function(a,T_term,X,Z_direct,Z_adaptive){
  N<-length(T_term);
  rho<-(1 - (1 - exp(-a * T_term))/(a * T_term))
  rhos<-matrix(rep(rho,Z_adaptive),nrow=N,ncol=Z_adaptive)
  dmX<-cbind(X[,1:Z_direct],X[,(Z_direct+1):(Z_adaptive+Z_direct)]*rhos)
  return(dmX)
}
