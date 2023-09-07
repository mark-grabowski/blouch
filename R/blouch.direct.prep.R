#' blouch.direct.prep - prep data to run in Stan's direct effect model
#'
#' @param trdata An object of the class treedata from function treeplyr
#' @param Y Vector containing name of column in treedata containing response variable
#' @param Y_error Vector containing name of column in treedata containing error of response variable
#' @param X Vector containing name(s) of column in treedata containing predictor variable(s)
#' @param X_error Vector containing name(s) of column in treedata containing error of predictor variable(s)
#'
#' @return dat - list file containing objecs setup for Blouch
#' @export
#'
blouch.direct.prep<-function(trdata,Y,Y_error,X,X_error){
  phy<-trdata$phy
  dat<-data.frame(trdata$dat)
  N<-length(trdata$phy$tip.label)
  ts<-ts_fxn(phy)
  ta<-ts[[1]]
  tij<-ts[[2]]
  datX<-as.matrix(dat[X])
  datXerror<-as.matrix(dat[X_error])
  Z_direct<-dim(X_error)[1]
  #print(paste(N,Z))
  #Direct effect model w/ Statistical Rethinking ME Correction
  dat<-list(N=N,Z_direct=Z,Y_obs=as.vector(t(dat[Y])),X_obs=matrix(datX,nrow=N,ncol=Z),
            Y_error=as.vector(t(dat[Y_error])),X_error=matrix(datXerror,nrow=N,ncol=Z),ta=ta,tij=tij)
  return(dat)
}
