#' sim.reg.direct.mlm.ve.data- Simulate Multilevel Multi-Optima Direct Effect Data - Varying Effects
#' @param phy An object of the class "phylo"
#' @param N Number of tips on tree
#' @param Z_direct Number of Direct Effect X traits
#' @param hl True half-life value
#' @param vy True Vy value
#' @param Sxx Instantaneous variance of the BM process
#' @param optima.bar True average optima
#' @param optima.sd True standard deviation of the optima
#' @param beta.bar True average beta
#' @param beta.sd True standard deviation of the betas
#' @param rho True corelation between optimas and betas
#' @param shifts Nodes for regime shifts
#' @return Merged phylogeny and data in treeplyr format
#' @export
#'
sim.reg.direct.mlm.ve.data<-function(phy,N,Z_direct,hl,vy,Sxx,optima.bar,optima.sd,beta.bar,beta.sd,rho,shifts){
  set.seed(10)
  N.regimes<-length(shifts)+1

  X<-phytools::fastBM(phy,a=0,sig2=Sxx,internal=FALSE) #Simulate X BM variable on tree, with scaling 10
  names(X)<-phy$tip.label
  phytools::phenogram(phy,X,spread.labels=TRUE,spread.cost=c(1,0)) #Plot X data

  trdata<-data.frame(phy$tip.label)
  trdata<-treeplyr::make.treedata(phy,trdata)
  trdata<-set.converge.regimes(trdata,shifts)

  #optima.sd<-2 #Standard deviation of the optima
  #beta.sd<-1 #Standard deviation of the slopes
  mlm.mu<-c(optima.bar,beta.bar) #Vector of means
  sigmas<-c(optima.sd,beta.sd) #Standard devations
  Rho<-matrix(c(1,rho,rho,1),nrow=2) #Correlation matrix
  Sigma<-diag(sigmas)%*%Rho%*%diag(sigmas) #Covariance matrix
  vary.effects<-MASS::mvrnorm(N.regimes,mlm.mu,Sigma) #Simulate regressions for each regimee
  optima.mlm<-vary.effects[,1] #Simulated optima/intercepts for regimes
  beta.mlm<-vary.effects[,2] #Simulate betas for each regime
  plot(optima.mlm,beta.mlm)
  vary.effects<-data.frame(Regimes=paste("OU",1:N.regimes,sep=""),Intercept=vary.effects[,1],Slope=vary.effects[,2])
  library(ellipse)
  for ( l in c(0.1,0.3,0.5,0.8,0.99)){
    lines(ellipse(Sigma,centre=mlm.mu,level=l))}

  reg_tips<-trdata$dat$regimes
  reg_tips<-as.numeric(as.factor(reg_tips))

  a<-log(2)/hl
  sigma2_y<-vy*(2*(log(2)/hl));

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
  lineages <- lapply(1:n, function(e) lineage.constructor(trdata$phy, e, anc_maps, regimes)) #Trace lineage from tips (n) to root and determine regimes of each node or branc

  optima.matrix<-weight.matrix(trdata$phy, a, lineages) #Slouch approach
  #return(list(optima.matrix,optima.mlm))
  mu<-matrix(NA,N,1)
  for(i in 1:(N)){
    mu[i]<-optima.matrix[i,]%*%optima.mlm+beta.mlm[reg_tips[i]]%*%X[i];
  }

  V<-calc_direct_V(phy,sigma2_y,a)
  Y<-MASS::mvrnorm(n=1,mu,V)

  X_error<-rep(0.01,N)
  Y_error<-rep(0.01,N)
  Y_with_error<-Y+rnorm(N,0,0.01)
  X_with_error<-X+rnorm(N,0,0.01)

  df <- data.frame( regime=reg_tips , Y=Y_with_error , X=X_with_error )

  ggplot2::ggplot(data=df,ggplot2::aes(y=Y_with_error,x=X_with_error,color=regime))+
    ggplot2::geom_point()

  trdata$dat<-cbind(trdata$dat,data.frame(cbind(reg_tips,Y_with_error,Y_error,X_with_error,X_error)))
  return(list(trdata,vary.effects))

}
#' sim.reg.adapt.mlm.ve.data- Simulate Multilevel Multi-Optima Adaptive Model Data - Varying Effects
#' @param phy An object of the class "phylo"
#' @param N Number of tips on tree
#' @param Z_adaptive Number of Adaptive X traits
#' @param hl True half-life value
#' @param vy True Vy value
#' @param Sxx Instantaneous variance of the BM process
#' @param optima.bar True average optima
#' @param optima.sd True standard deviation of the optima
#' @param beta.bar True average beta
#' @param beta.sd True standard deviation of the betas
#' @param rho True corelation between optimas and betas
#' @param shifts Nodes for regime shifts
#' @return Merged phylogeny and data in treeplyr format
#' @export
#'
sim.reg.adapt.mlm.ve.data<-function(phy,N,Z_adaptive,hl,vy,Sxx,optima.bar,optima.sd,beta.bar,beta.sd,rho,shifts){
  #Code to plot regimes
  set.seed(10)
  N.regimes<-length(shifts)+1
  trdata<-data.frame(phy$tip.label)
  trdata<-treeplyr::make.treedata(phy,trdata)
  trdata<-set.converge.regimes(trdata,shifts)

  shifts.total<-c(trdata$dat$regimes,trdata$phy$node.label)
  edge.regimes <- factor(shifts.total[trdata$phy$edge[,2]])

  #print(edge.regimes)
  reg.colors <- ggsci::pal_aaas("default", alpha = 0.7)(N.regimes)
  plot(trdata$phy,edge.color = reg.colors[edge.regimes], edge.width = 1, cex = 0.2)


  reg_tips<-trdata$dat$regimes
  reg_tips<-as.numeric(as.factor(reg_tips))

  print(reg.colors)
  plot(trdata$phy,edge.color = reg.colors[edge.regimes], edge.width = 1, cex = 0.2)


  #Calculate lineages
  n<-length(trdata$phy$tip.label)
  regimes_internal <-trdata$phy$node.label
  regimes_tip <- trdata$dat$regimes
  regimes <- concat.factor(regimes_tip, regimes_internal)
  anc_maps<-"regimes"
  lineages <- lapply(1:n, function(e) lineage.constructor(trdata$phy, e, anc_maps, regimes)) #Trace lineage from tips (n) to root and determine regimes of each node or branch

  ###############################################################################################################################################################################
  #MLM Code
  mlm.mu<-c(optima.bar,beta.bar) #Vector of means
  sigmas<-c(optima.sd,beta.sd) #Variance of intercepts and slopes

  Rho<-matrix(c(1,rho,rho,1),nrow=2) #Correlation matrix
  Sigma<-diag(sigmas)%*%Rho%*%diag(sigmas) #Covariance matrix - variance of intercepts and slopes, covariance of intercepts and slopes
  vary.effects<-MASS::mvrnorm(N.regimes,mlm.mu,Sigma) #Simulate regressions for each regimee
  optima.mlm<-vary.effects[,1] #Simulated optima/intercepts for regimes
  beta.mlm<-vary.effects[,2] #Simulate betas for each regime
  plot(optima.mlm,beta.mlm)
  vary.effects<-data.frame(Regimes=paste("OU",1:N.regimes,sep=""),Intercept=vary.effects[,1],Slope=vary.effects[,2])
  #Plot
  for ( l in c(0.1,0.3,0.5,0.8,0.99)){
    lines(ellipse::ellipse(Sigma,centre=mlm.mu,level=l))}

  ##################################################################################################################
  a<-log(2)/hl
  sigma2_y<-vy*(2*(log(2)/hl));

  vX0<-0
  vY0 <- 0

  X<-phytools::fastBM(phy,a=vX0,sig2=Sxx,internal=FALSE) #Simulate X BM variable on tree, with scaling 10
  names(X)<-phy$tip.label
  sigma2_x<-geiger::ratematrix(phy,X) #Calculate evolutionary v/cv matri
  phytools::phenogram(phy,X,spread.labels=TRUE,spread.cost=c(1,0)) #Plot X data

  optima_matrix<-weight.matrix(trdata$phy, a, lineages) #Slouch approach
  pred_X<-calc_adaptive_dmX(phy,a,X)

    mu<-matrix(NA,N,1)
  for(i in 1:N){
    #mu[i] = optima_matrix[i,]%*%optima+beta[reg_tips[i]]%*%pred_X[i]
    mu[i] = optima_matrix[i,]%*%optima.mlm+beta.mlm[reg_tips[i]]%*%pred_X[i]

  }

  n_reg<-length(unique(regimes))
  V<-calc_adaptive_V(phy,a, sigma2_y, beta.mlm,sigma2_x,Z_adaptive)
  Y<-MASS::mvrnorm(n=1,mu,V)
  ##################################################################################################################
  X_error<-rep(0.01,N)
  Y_error<-rep(0.01,N)
  Y_with_error<-Y+rnorm(N,0,0.01)
  X_with_error<-X+rnorm(N,0,0.01)

  trdata$dat<-cbind(trdata$dat,data.frame(cbind(reg_tips,Y_with_error,Y_error,X_with_error,X_error)))
  return(list(trdata,vary.effects))
}
