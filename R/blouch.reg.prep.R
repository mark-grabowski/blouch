#' blouch.reg.prep - setup dat file for Blouch's multi-optima model
#'
#' @param trdata An object of the class treedata from function treeplyr
#' @param Y Vector containing name of column in treedata containing response variable
#' @param Y_error Vector containing name of column in treedata containing error of response variable
#' @param reg.column Vector containing name of regime column in treedata$dat
#' @param anc_maps Vector containing name of regime type - at nodes "regimes" or SIMMAP
#'
#' @return dat - list file containing objecs setup for Blouch
#' @export
#'
blouch.reg.prep<-function(trdata,Y,Y_error,reg.column,anc_maps="regimes"){
  #Data should be send in treeplyr format
  #Only allows for regime shifts at nodes at present, not SIMMAP
  #Get Phylogeny info
  phy<-trdata$phy
  dat<-data.frame(trdata$dat)
  N<-length(trdata$phy$tip.label)

  n<-length(trdata$phy$tip.label)
  mrca1 <- ape::mrca(trdata$phy)
  times <- ape::node.depth.edgelength(trdata$phy)
  ta <- matrix(times[mrca1], nrow=n, dimnames = list(trdata$phy$tip.label, trdata$phy$tip.label))
  T.term <- times[1:n]
  tia <- times[1:n] - ta
  tja <- t(tia)
  tij <- tja + tia

  #Get internal regimes
  regimes_internal <-trdata$phy$node.label
  regimes_tip <- dat[reg.column][,1]

  regimes <- concat.factor(regimes_tip, regimes_internal)
  #anc_maps<-"regimes"
  lineages <- lapply(1:n, function(e) lineage.constructor(trdata$phy, e, anc_maps, regimes, ace)) #Trace lineage from tips (n) to root and determine regimes of each node or branch
  #return(lineages)
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
  #Y_obs<-trdata$dat$Y
  #Y_error<-trdata$dat$Y_error #Standard error not variance

  reg_tips<-dat[reg.column][,1]
  reg_tips<-as.numeric(as.factor(reg_tips))
  Dmat<-cophenetic(trdata$phy) #Time separating tips, same as tij matrix in Slouch/Blouch code

  ############################################################################################################
  #print(as.vector(t(dat[Y_error])))
  dat<-list(N=N,n_reg=length(unique(regimes)),max_node_num=max_node_num,Y_obs=as.vector(t(dat[Y])),Y_error=as.vector(t(dat[Y_error])),ta=ta,tij=tij,
            t_beginning=t_beginning,t_end=t_end,times=times,reg_match=reg_match,nodes=nodes,Dmat=Dmat)
  return(dat)
}
