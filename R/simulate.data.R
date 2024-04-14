##################################################################################################################
##################################################################################################################
#' sim.direct.data- Simulate Direct Effect Model Data
#' @param phy An object of the class "phylo"
#' @param N Number of tips on tree
#' @param Z Number of traits
#' @param hl True half-life value
#' @param vy True Vy value
#' @param Sxx Instantaneous variance of the BM process
#' @param optima True optima
#' @param beta True beta
#' @return Merged phylogeny and data in treeplyr format
#' @export
#'
sim.direct.data<-function(phy,N,Z,hl,vy,Sxx,optima,beta){
  set.seed(10)
  a<-log(2)/hl
  sigma2_y<-vy*(2*(log(2)/hl));
  vX0<-0
  vY0 <- 0
  ts<-ts_fxn(phy)
  ta<-ts[[1]]
  tij<-ts[[2]]
  V<-calc_direct_V(phy,sigma2_y,a)
  X<-phytools::fastBM(phy,a=vX0,sig2=Sxx,internal=FALSE) #Simulate X BM variable on tree, with scaling 10
  phytools::phenogram(phy,X,spread.labels=TRUE,spread.cost=c(1,0)) #Plot X data
  mu<-optima+X*beta #Simulate mu for Y
  #Simulate direct effect Y trait
  Y<-MASS::mvrnorm(n=1,mu,V)
  #plot
  df<-data.frame(Y=Y,X=X)
  names(df)<-c("Y","X")
  ggplot2::ggplot(data=df,ggplot2::aes(x=X,y=Y))+
    ggplot2::geom_point()
  summary(lm(Y~X,df))

  ########################################################################################################
  #Simulate errors - original Hansen setup
  X_error<-rep(0.01,N)
  Y_error<-rep(0.01,N)
  Y_with_error<-Y+rnorm(N,0,0.01)
  X_with_error<-X+rnorm(N,0,0.01)

  trait.data<-data.frame(cbind(Y_with_error,Y_error,X_with_error,X_error))
  trdata<-treeplyr::make.treedata(phy,trait.data)
  return(trdata)
}
##################################################################################################################
##################################################################################################################
#' sim.direct.multi.data- Simulate Direct Effect Model Data with Multiple X Traits
#' @param phy An object of the class "phylo"
#' @param N Number of tips on tree
#' @param Z Number of traits
#' @param hl True half-life value
#' @param vy True Vy value
#' @param Sxx Instantaneous variance of the BM process
#' @param optima True optima
#' @param beta True beta
#' @return Merged phylogeny and data in treeplyr format
#' @export
#'
sim.direct.multi.data<-function(phy,N,Z,hl,vy,Sxx,optima,beta){
  #Simulate multiple X traits - correlated
  set.seed(10)
  a<-log(2)/hl
  sigma2_y<-vy*(2*(log(2)/hl));
  vX0<-0
  vY0 <- 0
  ts<-ts_fxn(phy)
  ta<-ts[[1]]
  tij<-ts[[2]]
  V<-calc_direct_V(phy,sigma2_y,a)
  #vcv<-matrix(c(1,0.75,0.75,1),2,2) #Correlation between traits
  vcv<-matrix(c(1,0,0,1),2,2) #No correlation between traits
  Xs<-phytools::sim.corrs(phy,vcv) #Simulated correlated BM Xs
  #alpha<-2
  #beta<-c(0.35,0.1) #Slope
  mu<-optima+Xs%*%beta #Simulate mu for Y
  #Simulate direct effect Y trait
  Y<-mvrnorm(n=1,mu,V)
  #plot
  df<-data.frame(Y=Y,X=Xs)

  ggplot2::ggplot(data=df,ggplot2::aes(x=X.1,y=Y))+
    ggplot2::geom_point()
  summary(lm(Y~X.1,df))

  ggplot2::ggplot(data=df,ggplot2::aes(x=X.2,y=Y))+
    ggplot2::geom_point()
  summary(lm(Y~X.2,df))
  ########################################################################################################
  #Simulate errors with multiple traits - original Hansen setup
  X_error<-matrix(0.01,nrow=N,ncol=Z)
  Y_error<-rep(0.01,N)
  Y_with_error<-Y+rnorm(N,0,0.01)
  X_with_error<-apply(Xs,2,function(X){X+rnorm(N,0,0.01)})

  trait.data<-data.frame(cbind(Y_with_error,Y_error,X_with_error,X_error))
  trdata<-treeplyr::make.treedata(phy,trait.data)
  return(trdata)
}
##################################################################################################################
##################################################################################################################
#' sim.adaptive.data- Simulate Adaptive Data
#' @param phy An object of the class "phylo"
#' @param N Number of tips on tree
#' @param Z Number of traits
#' @param hl True half-life value
#' @param vy True Vy value
#' @param Sxx Instantaneous variance of the BM process
#' @param optima True optima
#' @param beta True beta
#' @return Merged phylogeny and data in treeplyr format
#' @export
#'
sim.adaptive.data<-function(phy,N,Z,hl,vy,Sxx,optima,beta){
  set.seed(10)
  ts<-ts_fxn(phy)
  ta<-ts[[1]]
  tij<-ts[[2]]
  T_term<-ts[[3]]
  tja<-ts[[4]]
  ########################################################################################################
  a<-log(2)/hl
  sigma2_y<-vy*(2*(log(2)/hl));
  vX0<-0
  vY0 <- 0
  X<-phytools::fastBM(phy,a=vX0,sig2=Sxx,internal=FALSE) #Simulate X BM variable on tree, with scaling 10
  phytools::phenogram(phy,X,spread.labels=TRUE,spread.cost=c(1,0)) #Plot X data
  sigma2_x<-geiger::ratematrix(phy,X) #Calculate evolutionary v/cv matrix
  #alpha<-2 #Intecept
  #beta<-0.25 #Slope
  dmX<-calc_adaptive_dmX(phy,a,X) #Calculate the design matrix
  mu<-optima+dmX%*%beta #Simulate mu for Y
  #V<-calc_adaptive_V(a, sigma2_y, ta,  tij,  tja,  T_term,  beta,  sigma2_x)
  V<-calc_adaptive_V(phy,a, sigma2_y, beta,sigma2_x,Z_adaptive)
  #Simulate direct effect Y trait
  Y<-MASS::mvrnorm(n=1,mu,V)
  #plot
  df<-data.frame(Y=Y,X=X)
  names(df)<-c("Y","X")

  ggplot2::ggplot(data=df,ggplot2::aes(x=X,y=Y))+
    ggplot2::geom_point()

  summary(lm(Y~X,df))

  ########################################################################################################
  #Simulate errors
  X_error<-rep(0.01,N)
  Y_error<-rep(0.01,N)
  Y_with_error<-Y+rnorm(N,0,0.01)
  X_with_error<-X+rnorm(N,0,0.01)

  #Make trdata file
  trait.data<-data.frame(cbind(Y_with_error,Y_error,X_with_error,X_error))
  trdata<-treeplyr::make.treedata(phy,trait.data)
  return(trdata)
}
##################################################################################################################
##################################################################################################################
#' sim.adaptive.multi.data- Simulate Adaptive Data with Multiple X Traits
#' @param phy An object of the class "phylo"
#' @param N Number of tips on tree
#' @param Z Number of traits
#' @param hl True half-life value
#' @param vy True Vy value
#' @param Sxx Instantaneous variance of the BM process
#' @param optima True optima
#' @param beta True beta
#' @return Merged phylogeny and data in treeplyr format
#' @export
#'
sim.adaptive.multi.data<-function(phy,N,Z,hl,vy,Sxx,optima,beta){
  set.seed(10)
  a<-log(2)/hl
  sigma2_y<-vy*(2*(log(2)/hl));
  vX0<-0
  vY0 <- 0
  vcv<-matrix(c(1,0,0,1),2,2) #No correlation between traits
  Xs<-phytools::sim.corrs(phy,vcv) #Simulated correlated BM Xs
  sigma2_x<-geiger::ratematrix(phy,Xs) #Calculate evolutionary v/cv matrix

  #beta<-c(0.35,0.1) #Slope
  dmX<-calc_dmX(a,T_term,Xs)
  mu<-alpha+dmX%*%beta #Simulate mu for Y

  V<-calc_adaptive_V(phy,a, sigma2_y, beta,sigma2_x,Z_adaptive)
  Y<-mvrnorm(n=1,mu,V)

  df<-data.frame(Y=Y,X=Xs)
  ggplot2::ggplot(data=df,ggplot2::aes(x=X.1,y=Y))+
    ggplot2::geom_point()

  summary(lm(Y~X.1,df))

  ggplot2::ggplot(data=df,ggplot2::aes(x=X.2,y=Y))+
    ggplot2::geom_point()

  summary(lm(Y~X.2,df))
  ########################################################################################################
  #Simulate errors with multiple traits - original Hansen setup
  X_error<-matrix(0.01,nrow=N,ncol=Z)
  Y_error<-rep(0.01,N)
  Y_with_error<-Y+rnorm(N,0,0.01)
  X_with_error<-apply(Xs,2,function(X){X+rnorm(N,0,0.01)})

  trdata$dat<-cbind(trdata$dat,data.frame(cbind(Y_with_error,Y_error,X_with_error,X_error)))
  return(trdata)
}
##################################################################################################################
##################################################################################################################
#' sim.adaptive.multi.data- Simulate Direct Effect + Adaptive Data with Multiple X Traits
#' @param phy An object of the class "phylo"
#' @param N Number of tips on tree
#' @param Z_direct Number of direct effect traits
#' @param Z_adaptive Number of adaptive traits
#' @param hl True half-life value
#' @param vy True Vy value
#' @param Sxx Instantaneous variance of the BM process
#' @param optima True optima
#' @param beta True beta
#' @return Merged phylogeny and data in treeplyr format
#' @export
#'
sim.direct.adaptive.data<-function(phy,N,Z_direct,Z_adaptive,hl,vy,Sxx,optima,beta){
  #Direct effect and Adaptive Model
  set.seed(10)
  #Setup parameters
  a<-log(2)/hl
  sigma2_y<-vy*(2*(log(2)/hl));

  vX0<-0
  vY0 <- 0

  ts<-ts_fxn(phy)
  ta<-ts[[1]]
  tij<-ts[[2]]
  T_term<-ts[[3]]
  tja<-ts[[4]]

  #Z_direct<-1
  #Z_adaptive<-1
  Z<-Z_direct+Z_adaptive

  #Xd<-rnorm(N,0,1)
  Xd<-phytools::fastBM(phy,a=vX0,sig2=Sxx,internal=FALSE) #Simulate X BM variable on tree, with scaling 10
  names(Xd)<-phy$tip.label
  phytools::phenogram(phy,Xd,spread.labels=TRUE,spread.cost=c(1,0)) #Plot X data
  Xa<-phytools::fastBM(phy,a=vX0,sig2=Sxx,internal=FALSE) #Simulate X BM variable on tree, with scaling 10
  names(Xa)<-phy$tip.label
  phytools::phenogram(phy,Xa,spread.labels=TRUE,spread.cost=c(1,0)) #Plot X data
  sigma2_x<-geiger::ratematrix(phy,Xa) #Calculate evolutionary v/cv matrix
  Xs<-cbind(Xd,Xa)

  #alpha<-2 #Intecept
  #beta<-c(0.35,0.25) #Slopes
  dmX<-calc_mixed_dmX(phy,a,Xs,Z_direct,Z_adaptive)
  mu<-optima+dmX%*%beta #Simulate mu for Y

  V<-calc_adaptive_V(phy,a, sigma2_y, beta[(Z_direct+1):(Z_adaptive+Z_direct)],sigma2_x,Z_adaptive)

  Y<-MASS::mvrnorm(n=1,mu,V)

  df<-data.frame(Y=Y,X=Xs)

  ggplot2::ggplot(data=df,ggplot2::aes(x=X.Xd,y=Y))+
    ggplot2::geom_point()

  summary(lm(Y~X.Xd,df))

  ggplot2::ggplot(data=df,ggplot2::aes(x=X.Xa,y=Y))+
    ggplot2::geom_point()

  summary(lm(Y~X.Xa,df))
  ############################################################################################################
  X_error<-matrix(0.01,nrow=N,ncol=Z)
  X_error<-data.frame(X_error)
  names(X_error)<-c("Xd_error","Xa_error")
  Y_error<-rep(0.01,N)
  Y_with_error<-Y+rnorm(N,0,0.01)
  X_with_error<-apply(Xs,2,function(X){X+rnorm(N,0,0.01)})
  ############################################################################################################
  #Make trdata file
  trait.data<-data.frame(cbind(Y_with_error,Y_error,X_with_error,X_error))
  trdata<-treeplyr::make.treedata(phy,trait.data)
  return(trdata)
}
##################################################################################################################
##################################################################################################################
#' sim.reg.data- Simulate Multi-Optima Data
#' @param phy An object of the class "phylo"
#' @param N Number of tips on tree
#' @param hl True half-life value
#' @param vy True Vy value
#' @param Sxx Instantaneous variance of the BM process
#' @param optima True optima
#' @param shifts Nodes for regime shifts
#' @return Merged phylogeny and data in treeplyr format
#' @export
#'
sim.reg.data<-function(phy,N,hl,vy,Sxx,optima,shits){
  set.seed(10)
  trdata<-data.frame(phy$tip.label)
  trdata<-treeplyr::make.treedata(phy,trdata)
  trdata<-set.converge.regimes(trdata,shifts)

  #Check if manual setting code worked
  shifts.total<-c(trdata$dat$regimes,trdata$phy$node.label)
  edge.regimes <- factor(shifts.total[trdata$phy$edge[,2]])
  print(edge.regimes)
  #Get ggplot colors used for plot to make on tree
  gg_color_hue <- function(n) {
    hues = seq(15, 375, length=n+1)
    hcl(h=hues, l=65, c=100)[1:n]
  }

  #reg.colors<-gg_color_hue(length(unique(trdata$dat$regimes)))
  reg.colors <- ggsci::pal_aaas("default", alpha = 0.7)(2)


  print(reg.colors)
  plot(trdata$phy,edge.color = reg.colors[edge.regimes], edge.width = 1, cex = 0.2)
  #Portrait 7X5

  #Simulate data
  n<-length(trdata$phy$tip.label)
  mrca1 <- ape::mrca(trdata$phy)
  times <- ape::node.depth.edgelength(trdata$phy)
  ta <- matrix(times[mrca1], nrow=n, dimnames = list(trdata$phy$tip.label, trdata$phy$tip.label))
  T.term <- times[1:n]
  tia <- times[1:n] - ta
  tja <- t(tia)
  tij <- tja + tia

  regimes_internal <-trdata$phy$node.label
  regimes_tip <- trdata$dat$regimes
  regimes <- concat.factor(regimes_tip, regimes_internal)
  anc_maps<-"regimes"
  lineages <- lapply(1:n, function(e) lineage.constructor(trdata$phy, e, anc_maps, regimes)) #Trace lineage from tips (n) to root and determine regimes of each node or branch

  ###################################################################################################################
  #Simulate Y based on V and incorporating regimes
  #Setup parameters

  #hl<-0.1 #0.1, 0.25, 0.75 - testing options
  a<-log(2)/hl
  #sigma2_y<-0.1
  #vy<-0.01
  sigma2_y<-vy*(2*(log(2)/hl))
  #alpha<-4 #Intecept
  #optima<-c(0.5,0.25) #Optima for two regimes

  dmX<-weight.matrix(trdata$phy, a, lineages) #Slouch approach
  mu<-dmX%*%optima #Simulate mu for Y
  V<-calc_direct_V(phy, sigma2_y, a)
  Y<-MASS::mvrnorm(n=1,mu,V)
  ########################################################################################################
  #Simulate errors with multiple traits - original Hansen setup
  Y_error<-rep(0.01,N)
  Y_with_error<-Y+rnorm(N,0,0.01)

  #Make trdata file
  trait.data<-data.frame(cbind(Y_with_error,Y_error))
  trdata$dat<-cbind(trdata$dat,data.frame(cbind(Y_with_error,Y_error)))
  return(trdata)
}
##################################################################################################################

##################################################################################################################
#' sim.reg.direct.data- Simulate Multi-Optima Direct Effect Model Data
#' @param phy An object of the class "phylo"
#' @param N Number of tips on tree
#' @param Z_direct Number of X traits
#' @param hl True half-life value
#' @param vy True Vy value
#' @param Sxx Instantaneous variance of the BM process
#' @param optima True optima
#' @param beta True beta
#' @param shifts Nodes for regime shifts
#' @return Merged phylogeny and data in treeplyr format
#' @export
#'
sim.reg.direct.data<-function(phy,N,Z_direct,hl,vy,Sxx,optima,beta,shifts){
  ########################################################################################################
  #Basic Setup - two regimes with direct only
  set.seed(10)
  trdata<-data.frame(phy$tip.label)
  trdata<-treeplyr::make.treedata(phy,trdata)
  trdata<-set.converge.regimes(trdata,shifts)

  #Check if manual setting code worked
  shifts.total<-c(trdata$dat$regimes,trdata$phy$node.label)
  edge.regimes <- factor(shifts.total[trdata$phy$edge[,2]])
  print(edge.regimes)

  reg.colors <- ggsci::pal_aaas("default", alpha = 0.7)(2)

  #print(reg.colors)
  plot(trdata$phy,edge.color = reg.colors[edge.regimes], edge.width = 1, cex = 0.2)

  #Phylogeny info
  n<-length(trdata$phy$tip.label)
  mrca1 <- ape::mrca(trdata$phy)
  times <- ape::node.depth.edgelength(trdata$phy)
  ta <- matrix(times[mrca1], nrow=n, dimnames = list(trdata$phy$tip.label, trdata$phy$tip.label))
  T.term <- times[1:n]
  tia <- times[1:n] - ta
  tja <- t(tia)
  tij <- tja + tia

  regimes_internal <-trdata$phy$node.label
  regimes_tip <- trdata$dat$regimes
  regimes <- concat.factor(regimes_tip, regimes_internal)
  anc_maps<-"regimes"
  lineages <- lapply(1:n, function(e) lineage.constructor(trdata$phy, e, anc_maps, regimes)) #Trace lineage from tips (n) to root and determine regimes of each node or branch

  #########################
  a<-log(2)/hl
  sigma2_y<-vy*(2*(log(2)/hl));

  vX0<-0
  vY0 <- 0
  Sxx<-10 #Look at effects

  ts<-ts_fxn(phy)
  ta<-ts[[1]]
  tij<-ts[[2]]
  T_term<-ts[[3]]
  tja<-ts[[4]]

  V<-calc_direct_V(phy,sigma2_y,a)
  X<-phytools::fastBM(phy,a=vX0,sig2=Sxx,internal=FALSE) #Simulate X BM variable on tree, with scaling 10
  #X<-rnorm(N,0,1)
  names(X)<-phy$tip.label
  phytools::phenogram(phy,X,spread.labels=TRUE,spread.cost=c(1,0)) #Plot X data
  dmX<-weight.matrix(trdata$phy, a, lineages) #Slouch approach
  dmX<-cbind(dmX,X)
  #beta<-c(2,1,0.25) #Two Optima/One Slope
  optima.beta<-c(optima,beta)
  mu<-dmX%*%optima.beta #Simulate mu for Y

  V<-calc_direct_V(phy,sigma2_y,a)
  Y<-MASS::mvrnorm(n=1,mu,V)

  ##################################################################################################################
  #Simulate errors - original Hansen setup
  X_error<-rep(0.01,N)
  Y_error<-rep(0.01,N)
  Y_with_error<-Y+rnorm(N,0,0.01)
  X_with_error<-X+rnorm(N,0,0.01)

  trdata$dat<-cbind(trdata$dat,data.frame(cbind(Y_with_error,Y_error,X_with_error,X_error)))
  return(trdata)
}
##################################################################################################################
##################################################################################################################
#' sim.reg.adapt.data- Simulate Multi-Optima Adaptive Model Data
#' @param phy An object of the class "phylo"
#' @param N Number of tips on tree
#' @param Z_adaptive Number of X traits
#' @param hl True half-life value
#' @param vy True Vy value
#' @param Sxx Instantaneous variance of the BM process
#' @param optima True optima
#' @param beta True beta
#' @param shifts Nodes for regime shifts
#' @return Merged phylogeny and data in treeplyr format
#' @export
#'
sim.reg.adapt.data<-function(phy,N,Z_adaptive,hl,vy,Sxx,optima,beta,shifts){
  ########################################################################################################
  #Basic Setup - two regimes with adaptive traits
  set.seed(10)
  trdata<-data.frame(phy$tip.label)
  trdata<-treeplyr::make.treedata(phy,trdata)
  trdata<-set.converge.regimes(trdata,shifts)

  #Check if manual setting code worked
  shifts.total<-c(trdata$dat$regimes,trdata$phy$node.label)
  edge.regimes <- factor(shifts.total[trdata$phy$edge[,2]])
  print(edge.regimes)

  reg.colors <- ggsci::pal_aaas("default", alpha = 0.7)(2)
  plot(trdata$phy,edge.color = reg.colors[edge.regimes], edge.width = 1, cex = 0.2)

  #Phylogeny info
  n<-length(trdata$phy$tip.label)
  mrca1 <- ape::mrca(trdata$phy)
  times <- ape::node.depth.edgelength(trdata$phy)
  ta <- matrix(times[mrca1], nrow=n, dimnames = list(trdata$phy$tip.label, trdata$phy$tip.label))
  T.term <- times[1:n]
  tia <- times[1:n] - ta
  tja <- t(tia)
  tij <- tja + tia

  regimes_internal <-trdata$phy$node.label
  regimes_tip <- trdata$dat$regimes
  regimes <- concat.factor(regimes_tip, regimes_internal)
  anc_maps<-"regimes"
  lineages <- lapply(1:n, function(e) lineage.constructor(trdata$phy, e, anc_maps, regimes)) #Trace lineage from tips (n) to root and determine regimes of each node or branch

  #########################
  a<-log(2)/hl
  sigma2_y<-vy*(2*(log(2)/hl));
  vX0<-0
  vY0 <- 0

  ts<-ts_fxn(phy)
  ta<-ts[[1]]
  tij<-ts[[2]]
  T_term<-ts[[3]]
  tja<-ts[[4]]

  X<-phytools::fastBM(phy,a=vX0,sig2=Sxx,internal=FALSE) #Simulate X BM variable on tree, with scaling 10
  names(X)<-phy$tip.label
  phytools::phenogram(phy,X,spread.labels=TRUE,spread.cost=c(1,0)) #Plot X data
  sigma2_x<-geiger::ratematrix(phy,X) #Calculate evolutionary v/cv matrix
  dmX<-weight.matrix(trdata$phy, a, lineages) #Slouch approach
  dmX<-cbind(dmX,calc_adaptive_dmX(phy,a,X))
  optima.beta<-c(optima,beta)
  #beta<-c(2,1,0.25) #Two Optima/Two Slopes
  mu<-dmX%*%optima.beta #Simulate mu for Y

  V<-calc_adaptive_V(phy,a, sigma2_y, beta,sigma2_x,Z_adaptive)
  Y<-MASS::mvrnorm(n=1,mu,V)
#######################################################################################################
  X_error<-rep(0.01,N)
  Y_error<-rep(0.01,N)
  Y_with_error<-Y+rnorm(N,0,0.01)
  X_with_error<-X+rnorm(N,0,0.01)

  #Make trdata file
  trdata$dat<-cbind(trdata$dat,data.frame(cbind(Y_with_error,Y_error,X_with_error,X_error)))
  return(trdata)
}
#' sim.reg.direct.ve.data- Simulate Multi-Optima Direct Effect Model Data - Varying Effects
#' @param phy An object of the class "phylo"
#' @param N Number of tips on tree
#' @param Z_direct Number of X traits
#' @param hl True half-life value
#' @param vy True Vy value
#' @param Sxx Instantaneous variance of the BM process
#' @param optima True optima
#' @param beta True beta
#' @param shifts Nodes for regime shifts
#' @return Merged phylogeny and data in treeplyr format
#' @export
#'
sim.reg.direct.ve.data<-function(phy,N,Z_direct,hl,vy,Sxx,optima,beta,shifts){
  set.seed(10)
  trdata<-data.frame(phy$tip.label)
  trdata<-treeplyr::make.treedata(phy,trdata)
  trdata<-set.converge.regimes(trdata,shifts)

  #Check if manual setting code worked
  shifts.total<-c(trdata$dat$regimes,trdata$phy$node.label)
  edge.regimes <- factor(shifts.total[trdata$phy$edge[,2]])
  print(edge.regimes)

  reg.colors <- ggsci::pal_aaas("default", alpha = 0.7)(2)
  plot(trdata$phy,edge.color = reg.colors[edge.regimes], edge.width = 1, cex = 0.2)

  reg_tips<-trdata$dat$regimes
  reg_tips<-as.numeric(as.factor(reg_tips))

  print(reg.colors)
  plot(trdata$phy,edge.color = reg.colors[edge.regimes], edge.width = 1, cex = 0.2)

  #Phylogeny info
  n<-length(trdata$phy$tip.label)
  mrca1 <- ape::mrca(trdata$phy)
  times <- ape::node.depth.edgelength(trdata$phy)
  ta <- matrix(times[mrca1], nrow=n, dimnames = list(trdata$phy$tip.label, trdata$phy$tip.label))
  T.term <- times[1:n]
  tia <- times[1:n] - ta
  tja <- t(tia)
  tij <- tja + tia

  regimes_internal <-trdata$phy$node.label
  regimes_tip <- trdata$dat$regimes
  regimes <- concat.factor(regimes_tip, regimes_internal)
  anc_maps<-"regimes"
  lineages <- lapply(1:n, function(e) lineage.constructor(trdata$phy, e, anc_maps, regimes)) #Trace lineage from tips (n) to root and determine regimes of each node or branch

  #########################
  #hl<-0.25 #0.1, 0.25, 0.75 - testing options
  a<-log(2)/hl
  #vy<-0.01 #0.25,0.5 - testing options
  sigma2_y<-vy*(2*(log(2)/hl));

  vX0<-0
  vY0 <- 0
  #Sxx<-10 #Look at effects

  ts<-ts_fxn(phy)
  ta<-ts[[1]]
  tij<-ts[[2]]
  T_term<-ts[[3]]
  tja<-ts[[4]]


  X<-phytools::fastBM(phy,a=vX0,sig2=Sxx,internal=FALSE) #Simulate X BM variable on tree, with scaling 10
  names(X)<-phy$tip.label
  phytools::phenogram(phy,X,spread.labels=TRUE,spread.cost=c(1,0)) #Plot X data

  dmX<-weight.matrix(trdata$phy, a, lineages) #Slouch approach
  #optima<-c(2,1)
  #beta<-c(0.25,0.15) #Two Optima/Two Slopes
  mu<-matrix(NA,N,1)
  for(i in 1:N){
    mu[i]<-dmX[i,]%*%optima+beta[reg_tips[i]]%*%X[i];
  }

  V<-calc_direct_V(phy,sigma2_y,a)
  Y<-MASS::mvrnorm(n=1,mu,V)

  X_error<-rep(0.01,N)
  Y_error<-rep(0.01,N)
  Y_with_error<-Y+rnorm(N,0,0.01)
  X_with_error<-X+rnorm(N,0,0.01)

  trdata$dat<-cbind(trdata$dat,data.frame(cbind(Y_with_error,Y_error,X_with_error,X_error)))
  return(trdata)

}
#' sim.reg.adapt.ve.data- Simulate Multi-Optima Adaptive Model Data - Varying Effects
#' @param phy An object of the class "phylo"
#' @param N Number of tips on tree
#' @param Z_adaptive Number of X traits
#' @param hl True half-life value
#' @param vy True Vy value
#' @param Sxx Instantaneous variance of the BM process
#' @param optima True optima
#' @param beta True beta
#' @param shifts Nodes for regime shifts
#' @return Merged phylogeny and data in treeplyr format
#' @export
#'
sim.reg.adapt.ve.data<-function(phy,N,Z_adaptive,hl,vy,Sxx,optima,beta,shifts){
  set.seed(10)
  trdata<-data.frame(phy$tip.label)
  trdata<-treeplyr::make.treedata(phy,trdata)
  trdata<-set.converge.regimes(trdata,shifts)
  #Check if manual setting code worked
  shifts.total<-c(trdata$dat$regimes,trdata$phy$node.label)
  edge.regimes <- factor(shifts.total[trdata$phy$edge[,2]])
  print(edge.regimes)
  reg.colors <- ggsci::pal_aaas("default", alpha = 0.7)(2)
  plot(trdata$phy,edge.color = reg.colors[edge.regimes], edge.width = 1, cex = 0.2)

  reg_tips<-trdata$dat$regimes
  reg_tips<-as.numeric(as.factor(reg_tips))

  print(reg.colors)
  plot(trdata$phy,edge.color = reg.colors[edge.regimes], edge.width = 1, cex = 0.2)

  #Phylogeny info
  n<-length(trdata$phy$tip.label)
  mrca1 <- ape::mrca(trdata$phy)
  times <- ape::node.depth.edgelength(trdata$phy)
  ta <- matrix(times[mrca1], nrow=n, dimnames = list(trdata$phy$tip.label, trdata$phy$tip.label))
  T.term <- times[1:n]
  tia <- times[1:n] - ta
  tja <- t(tia)
  tij <- tja + tia

  regimes_internal <-trdata$phy$node.label
  regimes_tip <- trdata$dat$regimes
  regimes <- concat.factor(regimes_tip, regimes_internal)
  anc_maps<-"regimes"
  lineages <- lapply(1:n, function(e) lineage.constructor(trdata$phy, e, anc_maps, regimes)) #Trace lineage from tips (n) to root and determine regimes of each node or branch

  ###############################################################################################################################################################################
  a<-log(2)/hl
  sigma2_y<-vy*(2*(log(2)/hl));

  vX0<-0
  vY0 <- 0

  ts<-ts_fxn(phy)
  ta<-ts[[1]]
  tij<-ts[[2]]
  T_term<-ts[[3]]
  tja<-ts[[4]]

  X<-phytools::fastBM(phy,a=vX0,sig2=Sxx,internal=FALSE) #Simulate X BM variable on tree, with scaling 10
  names(X)<-phy$tip.label
  sigma2_x<-geiger::ratematrix(phy,X) #Calculate evolutionary v/cv matri
  phytools::phenogram(phy,X,spread.labels=TRUE,spread.cost=c(1,0)) #Plot X data

  optima_matrix<-weight.matrix(trdata$phy, a, lineages) #Slouch approach
  pred_X<-calc_adaptive_dmX(phy,a,X)
  #dmX<-weight.matrix(trdata$phy, a, lineages) #Slouch approach
  #dmX<-cbind(dmX,calc_adaptive_dmX(phy,a,X))
  #optima<-c(2,1)
  #beta<-c(0.25,0.15) #Two Optima/Two Slopes
  mu<-matrix(NA,N,1)
  for(i in 1:N){
    mu[i] = optima_matrix[i,]%*%optima+beta[reg_tips[i]]%*%pred_X[i]
  }

  n_reg<-length(unique(regimes))
  V<-calc_adaptive_V(phy,a, sigma2_y, beta,sigma2_x,Z_adaptive)
  Y<-MASS::mvrnorm(n=1,mu,V)
  ##################################################################################################################
  X_error<-rep(0.01,N)
  Y_error<-rep(0.01,N)
  Y_with_error<-Y+rnorm(N,0,0.01)
  X_with_error<-X+rnorm(N,0,0.01)

  trdata$dat<-cbind(trdata$dat,data.frame(cbind(Y_with_error,Y_error,X_with_error,X_error)))
  return(trdata)
}
#' sim.reg.direct.adapt.data- Simulate Multi-Optima Direct EFfect + Adaptive Model Data
#' @param phy An object of the class "phylo"
#' @param N Number of tips on tree
#' @param Z_direct Number of Direct Effect X traits
#' @param Z_adaptive Number of Adaptive X traits
#' @param hl True half-life value
#' @param vy True Vy value
#' @param Sxx Instantaneous variance of the BM process
#' @param optima True optima
#' @param beta True beta
#' @param shifts Nodes for regime shifts
#' @return Merged phylogeny and data in treeplyr format
#' @export
#'
sim.reg.direct.adapt.data<-function(phy,N,Z_direct,Z_adaptive,hl,vy,Sxx,optima,beta,shifts){
  set.seed(10)
  trdata<-data.frame(phy$tip.label)
  trdata<-treeplyr::make.treedata(phy,trdata)
  trdata<-set.converge.regimes(trdata,shifts)
  #Check if manual setting code worked
  shifts.total<-c(trdata$dat$regimes,trdata$phy$node.label)
  edge.regimes <- factor(shifts.total[trdata$phy$edge[,2]])
  print(edge.regimes)

  reg.colors <- ggsci::pal_aaas("default", alpha = 0.7)(2)
  plot(trdata$phy,edge.color = reg.colors[edge.regimes], edge.width = 1, cex = 0.2)

  print(reg.colors)
  plot(trdata$phy,edge.color = reg.colors[edge.regimes], edge.width = 1, cex = 0.2)

  #Phylogeny info
  n<-length(trdata$phy$tip.label)
  mrca1 <- ape::mrca(trdata$phy)
  times <- ape::node.depth.edgelength(trdata$phy)
  ta <- matrix(times[mrca1], nrow=n, dimnames = list(trdata$phy$tip.label, trdata$phy$tip.label))
  T.term <- times[1:n]
  tia <- times[1:n] - ta
  tja <- t(tia)
  tij <- tja + tia

  regimes_internal <-trdata$phy$node.label
  regimes_tip <- trdata$dat$regimes
  regimes <- concat.factor(regimes_tip, regimes_internal)
  anc_maps<-"regimes"
  lineages <- lapply(1:n, function(e) lineage.constructor(trdata$phy, e, anc_maps, regimes)) #Trace lineage from tips (n) to root and determine regimes of each node or branch

  #########################
  a<-log(2)/hl
  sigma2_y<-vy*(2*(log(2)/hl));

  vX0<-0
  vY0 <- 0

  ts<-ts_fxn(phy)
  ta<-ts[[1]]
  tij<-ts[[2]]
  T_term<-ts[[3]]
  tja<-ts[[4]]
  Z<-Z_direct+Z_adaptive

  #Xd<-rnorm(N,0,1)
  Xd<-phytools::fastBM(phy,a=vX0,sig2=Sxx,internal=FALSE) #Simulate X BM variable on tree, with scaling 10
  names(Xd)<-phy$tip.label
  phytools::phenogram(phy,Xd,spread.labels=TRUE,spread.cost=c(1,0)) #Plot X data
  Xa<-phytools::fastBM(phy,a=vX0,sig2=Sxx,internal=FALSE) #Simulate X BM variable on tree, with scaling 10
  names(Xa)<-phy$tip.label
  phytools::phenogram(phy,Xa,spread.labels=TRUE,spread.cost=c(1,0)) #Plot X data
  sigma2_x<-geiger::ratematrix(phy,Xa) #Calculate evolutionary v/cv matrix
  Xs<-cbind(Xd,Xa)

  dmX<-weight.matrix(trdata$phy, a, lineages) #Slouch approach
  dmX<-cbind(dmX,calc_mixed_dmX(phy,a,Xs,Z_direct,Z_adaptive))
  #beta<-c(2,1,0.35,0.25) #Two Optima/Two Slopes
  optima.beta<-c(optima,beta)
  mu<-dmX%*%optima.beta #Simulate mu for Y

  V<-calc_adaptive_V(phy,a, sigma2_y, beta[Z_direct+1],  sigma2_x, Z_adaptive)
  Y<-MASS::mvrnorm(n=1,mu,V)

  ########################################################################################################
  #Simulate errors - for use with blouchOU_reg_direct_adaptive_ME
  X_error<-matrix(0.01,nrow=N,ncol=Z)
  Y_error<-rep(0.01,N)
  Y_with_error<-Y+rnorm(N,0,0.01)
  X_with_error<-apply(Xs,2,function(X){X+rnorm(N,0,0.01)})
  trdata$dat<-cbind(trdata$dat,data.frame(cbind(Y_with_error,Y_error,X_with_error,X_error)))
  names(trdata$dat)[6:7]<-c("Xd_error","Xa_error")
  return(trdata)

}
#' sim.reg.direct.adapt.data- Simulate Multi-Optima Direct EFfect + Adaptive Model Data - Varying Effects
#' @param phy An object of the class "phylo"
#' @param N Number of tips on tree
#' @param Z_direct Number of Direct Effect X traits
#' @param Z_adaptive Number of Adaptive X traits
#' @param hl True half-life value
#' @param vy True Vy value
#' @param Sxx Instantaneous variance of the BM process
#' @param optima True optima
#' @param beta True beta
#' @param shifts Nodes for regime shifts
#' @return Merged phylogeny and data in treeplyr format
#' @export
#'
sim.reg.direct.adapt.ve.data<-function(phy,N,Z_direct,Z_adaptive,hl,vy,Sxx,optima,beta,shifts){
  set.seed(10)
  #Four regimes with 1 direct and 1 adaptive trait and multiple slopes per optima but single alpha parameter
  trdata<-data.frame(phy$tip.label)
  trdata<-treeplyr::make.treedata(phy,trdata)
  trdata<-set.converge.regimes(trdata,shifts)

  #Check if manual setting code worked
  shifts.total<-c(trdata$dat$regimes,trdata$phy$node.label)
  edge.regimes <- factor(shifts.total[trdata$phy$edge[,2]])
  print(edge.regimes)

  reg.colors <- ggsci::pal_aaas("default", alpha = 0.7)(2)
  plot(trdata$phy,edge.color = reg.colors[edge.regimes], edge.width = 1, cex = 0.2)

  reg_tips<-trdata$dat$regimes
  reg_tips<-as.numeric(as.factor(reg_tips))

  print(reg.colors)
  plot(trdata$phy,edge.color = reg.colors[edge.regimes], edge.width = 1, cex = 0.2)

  #Phylogeny info
  n<-length(trdata$phy$tip.label)
  mrca1 <- ape::mrca(trdata$phy)
  times <- ape::node.depth.edgelength(trdata$phy)
  ta <- matrix(times[mrca1], nrow=n, dimnames = list(trdata$phy$tip.label, trdata$phy$tip.label))
  T.term <- times[1:n]
  tia <- times[1:n] - ta
  tja <- t(tia)
  tij <- tja + tia
  Z<-Z_direct+Z_adaptive

  regimes_internal <-trdata$phy$node.label
  regimes_tip <- trdata$dat$regimes
  regimes <- concat.factor(regimes_tip, regimes_internal)
  anc_maps<-"regimes"
  lineages <- lapply(1:n, function(e) lineage.constructor(trdata$phy, e, anc_maps, regimes)) #Trace lineage from tips (n) to root and determine regimes of each node or branch

  ###############################################################################################################################################################################
  a<-log(2)/hl
  sigma2_y<-vy*(2*(log(2)/hl));

  vX0<-0
  vY0 <- 0

  ts<-ts_fxn(phy)
  ta<-ts[[1]]
  tij<-ts[[2]]
  T_term<-ts[[3]]
  tja<-ts[[4]]

  Xa<-phytools::fastBM(phy,a=vX0,sig2=Sxx,internal=FALSE) #Simulate X BM variable on tree, with scaling 10
  Xd<-phytools::fastBM(phy,a=vX0,sig2=Sxx,internal=FALSE) #Simulate X BM variable on tree, with scaling 10
  names(Xa)<-phy$tip.label
  names(Xd)<-phy$tip.label
  phytools::phenogram(phy,Xd,spread.labels=TRUE,spread.cost=c(1,0)) #Plot X data
  phytools::phenogram(phy,Xa,spread.labels=TRUE,spread.cost=c(1,0)) #Plot X data

  Xs<-cbind(Xd,Xa)
  sigma2_x<-geiger::ratematrix(phy,Xa) #Calculate evolutionary v/cv matrix

  optima_matrix<-weight.matrix(trdata$phy, a, lineages) #Slouch approach
  pred_X<-calc_mixed_dmX(phy,a,Xs,Z_direct,Z_adaptive)

  #optima<-c(2,1)
  #beta<-c(0.25,0.15,0.35,0.1) #Two Optima/Two Slopes
  #beta<-data.frame(matrix(c(0.25,0.15,0.35,0.1),ncol=2,nrow=2)) #Two traits on columns, two regimes on vertical

  mu<-matrix(NA,N,1)
  for(i in 1:N){
    #  mu[i] = optima_matrix[i,]%*%optima+beta[reg_tips[i]]%*%pred_X[i]
    mu[i] = optima_matrix[i,]%*%optima+pred_X[i,]%*%t(beta[reg_tips[i],])

  }

  n_reg<-length(unique(regimes))
  V<-calc_adaptive_V(phy,a, sigma2_y, beta[Z_direct+1],  sigma2_x, Z_adaptive)
  Y<-MASS::mvrnorm(n=1,mu,V)

  ##################################################################################################################
  X_error<-matrix(0.01,nrow=N,ncol=2)
  Y_error<-rep(0.01,N)
  Y_with_error<-Y+rnorm(N,0,0.01)
  X_with_error<-Xs+rnorm(N,0,0.01)
  trdata$dat<-cbind(trdata$dat,data.frame(cbind(Y_with_error,Y_error,X_with_error,X_error)))
  names(trdata$dat)[6:7]<-c("Xd_error","Xa_error")
  return(trdata)
}
