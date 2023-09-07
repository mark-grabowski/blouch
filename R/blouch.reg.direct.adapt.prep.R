#' blouch.reg.direct.adapt.prep - function to setup dat file for Blouch's multi-optima direct effect adaptive model
#' @param trdata An object of the class treedata from function treeplyr
#' @param Y Vector containing name of column in treedata containing response variable
#' @param Y_error Vector containing name of column in treedata containing error of response variable
#' @param X Vector containing name(s) of column in treedata containing predictor variable(s)
#' @param X_error Vector containing name(s) of column in treedata containing error of predictor variable(s)
#' @param Z_direct Vector containing number of direct effect predictor traits
#' @param Z_adaptive Vector containing number of adaptive predictor traits
#' @param reg.column Vector containing name of regime column in treedata$dat
#'
#' @return dat - list file containing objecs setup for Blouch
#' @export
#'
blouch.reg.direct.adapt.prep<-function(trdata,Y,Y_error,X,X_error,Z_direct,Z_adaptive,reg.column){
  #Data should be send in treeplyr format
  #Only allows for regime shifts at nodes at present, not SIMMAP
  #Get Phylogeny info
  anc_maps="regimes"
  phy<-trdata$phy
  dat<-data.frame(trdata$dat)
  #return(dat)
  N<-length(trdata$phy$tip.label)

  n<-length(trdata$phy$tip.label)
  mrca1 <- ape::mrca(trdata$phy)
  times <- ape::node.depth.edgelength(trdata$phy)
  ta <- matrix(times[mrca1], nrow=n, dimnames = list(trdata$phy$tip.label, trdata$phy$tip.label))
  T_term <- times[1:n]
  tia <- times[1:n] - ta
  tja <- t(tia)
  tij <- tja + tia

  #Get internal regimes
  regimes_internal <-trdata$phy$node.label
  regimes_tip <- dat[reg.column][,1]
  regimes <- concat.factor(regimes_tip, regimes_internal)
  #anc_maps<-"regimes"
  lineages <- lapply(1:n, function(e) lineage.constructor(trdata$phy, e, anc_maps, regimes, ace)) #Trace lineage from tips (n) to root and determine regimes of each node or branch
  ############################################################################################################
  nodes<-NULL
  store<-NULL
  reg_num_lineage<-NULL
  for(i in 1:length(lineages)){
    store<-c(store,length(lineage.nodes(trdata$phy,i))) #Calculate max node height
    reg_num_lineage<-c(reg_num_lineage,length(unique(lineages[[i]]$lineage_regimes)))
    nodes<-c(nodes,length(lineages[[i]]$nodes))
  }
  max_node_num<-max(store)
  times<-matrix(0,length(lineages),max_node_num)
  t_end<-matrix(0,length(lineages),max_node_num)
  t_beginning<-matrix(0,length(lineages),max_node_num)
  reg_match<-data.frame(matrix(0,length(lineages),max_node_num))

  for(i in 1:length(lineages)){
    times[i,1:length(lineages[[i]]$times)]<-lineages[[i]]$times
    t_end[i,1:length(lineages[[i]]$t_end)]<-lineages[[i]]$t_end
    t_beginning[i,1:length(lineages[[i]]$t_beginning)]<-lineages[[i]]$t_beginning
    reg_match[i,1:length(lineages[[i]]$lineage_regimes)]<-rev(as.numeric(lineages[[i]]$lineage_regimes))
  }
  ############################################################################################################
  Z_X_error<-length(X_error) #Number of X traits with error

  reg_tips<-dat[reg.column][,1]
  reg_tips<-as.numeric(as.factor(reg_tips))
  Dmat<-cophenetic(trdata$phy) #Time separating tips, same as tij matrix in Slouch/Blouch code
  ############################################################################################################

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

  dat<-list(N=N,n_reg=length(unique(regimes)),Z_direct=Z_direct,Z_adaptive=Z_adaptive,Z_X_error=Z_X_error,max_node_num=max_node_num,
            Y_obs=as.vector(t(dat[Y])),X_obs=data.matrix(dat[X]),Y_error=as.vector(t(dat[Y_error])),
            X_error=data.matrix(dat[X_error]),sigma2_x=sigma2_x,ta=ta,tij=tij,
            tja=tja,T_term=T_term,t_beginning=t_beginning,
            t_end=t_end,times=times,reg_match=reg_match,nodes=nodes,reg_tips=reg_tips)


  return(dat)
}
