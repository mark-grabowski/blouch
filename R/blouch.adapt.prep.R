#' blouch.adapt.prep - setup dat file to run in Blouch's adaptive model
#'
#' @param trdata An object of the class treedata from function treeplyr
#' @param Y Vector containing name of column in treedata containing response variable
#' @param Y_error Vector containing name of column in treedata containing error of response variable
#' @param X Vector containing name(s) of column in treedata containing predictor variable(s)
#' @param X_error Vector containing name(s) of column in treedata containing error of predictor variable(s)
#' @param Z_adaptive Number of adaptive traits
#'
#' @return dat - list file containing objecs setup for Blouch
#' @export
#'
blouch.adapt.prep<-function(trdata,Y,Y_error,X,X_error,Z_adaptive){
  phy<-trdata$phy
  dat<-data.frame(trdata$dat)
  N<-length(trdata$phy$tip.label)
  ts<-ts_fxn(phy)
  ta<-ts[[1]]
  tij<-ts[[2]]
  T_term<-ts[[3]]
  tja<-ts[[4]]

  datX<-as.matrix(dat[X])
  datXerror<-as.matrix(dat[X_error])

  if(is.na(Y_error)!=TRUE){mv.response<-dat[Y_error]}
  else{mv.response<-as.vector(rep(0,N))}
  if(is.na(X_error)!=TRUE){mv.pred<-matrix(datXerror,nrow=N,ncol=Z_adaptive)}
  else{mv.pred<-matrix(0,nrow=N,ncol=Z_adaptive)}

  test<-sigma.X.estimate(phy,ta, predictor = datX, mv.predictor = mv.pred)

  brownian_mean<-test[1]
  sigma_squared_x<-test[2]

  #Direct effect model w/ Statistical Rethinking ME Correction
  dat<-list(N=N,Z_adapt=Z_adaptive,Y_obs=as.vector(t(dat[Y])),X_obs=matrix(datX,nrow=N,ncol=Z_adaptive),
            Y_error=as.vector(t(dat[Y_error])),X_error=mv.pred,
            ta=ta,tij=tij,tja=tja,T_term=T_term,sigma2_x=sigma_squared_x)
  return(dat)
}
