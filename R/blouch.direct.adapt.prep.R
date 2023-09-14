#' blouch.direct.adapt.prep - setup dat file to run in Blouch's direct effect and adpative model
#'
#' @param trdata An object of the class treedata from function treeplyr
#' @param Y Vector containing name of column in treedata containing response variable
#' @param Y_error Vector containing name of column in treedata containing error of response variable
#' @param X Vector containing name(s) of column in treedata containing predictor variable(s)
#' @param X_error Vector containing name(s) of column in treedata containing error of predictor variable(s)
#' @param Z_direct Vector containing number of direct effect predictor traits
#' @param Z_adaptive Vector containing number of adaptive predictor traits
#'
#' @return dat - list file containing objecs setup for Blouch
#' @export
#'
blouch.direct.adapt.prep<-function(trdata,Y,Y_error,X,X_error,Z_direct,Z_adaptive){
  phy<-trdata$phy
  dat<-data.frame(trdata$dat)
  N<-length(trdata$phy$tip.label)
  ts<-ts_fxn(phy)
  ta<-ts[[1]]
  tij<-ts[[2]]
  T_term<-ts[[3]]
  tja<-ts[[4]]

  datX<-as.matrix(dat[X])
  #return(Z_adapt)
  datX.direct<-datX[,1:Z_direct]
  datX.adapt<-datX[,(Z_direct+1):(Z_adaptive+Z_direct)]

  datXerror<-as.matrix(dat[X_error])

  datXerror.direct<-datXerror[,1:Z_direct]
  datXerror.adapt<-datXerror[,(Z_direct+1):(Z_adaptive+Z_direct)]
  #return(datXerror.direct)
  if(is.na(Y_error)!=TRUE){mv.response<-dat[Y_error]}
  else{mv.response<-as.vector(rep(0,N))}

  if(any(is.na(datXerror.direct)!=TRUE)){mv.pred.direct<-matrix(datXerror.direct,nrow=N,ncol=Z_direct)}
  else{mv.pred.direct<-matrix(0,nrow=N,ncol=Z_direct)}

  if(any(is.na(datXerror.adapt)!=TRUE)){mv.pred.adapt<-matrix(datXerror.adapt,nrow=N,ncol=Z_adaptive)}
  else{mv.pred.adapt<-matrix(0,nrow=N,ncol=Z_adaptive)}

  test<-sigma.X.estimate(phy,ta, predictor = datX.adapt, mv.predictor = mv.pred.adapt)

  brownian_mean<-test[1]
  sigma_squared_x<-test[2]

  Z<-Z_direct+Z_adaptive
  Z_X_error<-Z_direct+Z_adaptive

  dat<-list(N=N,Z_direct=Z_direct,Z_adaptive=Z_adaptive,Z_X_error=Z_X_error,
            Y_obs=as.vector(t(dat[Y])),X_obs=matrix(datX,nrow=N,ncol=Z),
            Y_error=as.vector(t(dat[Y_error])),X_error=matrix(datXerror,nrow=N,ncol=Z),
            ta=ta,tij=tij,tja=tja,T_term=T_term,sigma2_x=sigma_squared_x)
  return(dat)
}
