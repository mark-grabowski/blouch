
#' set.converge.regimes - function to assign regimes on a phylogeny
#'
#' @param trdata treeplyr format file
#' @param regimes node numbers for where the regimes shift
#'
#' @return treeplyr format file with new column in dat for tip regimes and internal regime assignments in node.labels
#' @export
#'
set.converge.regimes<-function(trdata,regimes){
  #02.19.23 - allows convergence in regimes
  #returns a new column of the dataset - regimes
  #also returns internal node assignment
  #Make sure you send a merged trdata file from treeplyr
  getDescendants<-function(tree,node,curr=NULL){
    if(is.null(curr)) curr<-vector()
    daughters<-tree$edge[which(tree$edge[,1]==node),2]
    curr<-c(curr,daughters)
    w<-which(daughters>=length(tree$tip))
    if(length(w)>0) for(i in 1:length(w))
      curr<-getDescendants(tree,daughters[w[i]],curr)
    return(curr)
  }

  n.tips<-length(trdata$phy$tip.label)
  num.internal.nodes<-n.tips+(1:trdata$phy$Nnode)
  n.internal.nodes<-length(num.internal.nodes)

  #First set all exant tips to OU1
  trdata$dat$regimes<-"OU1"
  #Assign all internal nodes to OU1
  internal.nodes.regimes<-rep("OU1",n.internal.nodes)

  for(i in 1:(length(regimes))){
    print(i)
    for(j in 1:length(regimes[[i]])){
      if(regimes[[i]][[j]]<=n.tips){ #shift on single branch
        trdata$dat$regimes[regimes[[i]][[j]]]<-paste("OU",sep="",i+1 )#Assign external nodes/tips to OUi
      }
      for(j in 1:length(regimes[[i]])){
        rep.regimes<-regimes[[i]][[j]]
        saved.decendants<-getDescendants(trdata$phy,node=rep.regimes)
        external.nodes<-saved.decendants[saved.decendants<=n.tips]
        internal.nodes<-saved.decendants[saved.decendants>n.tips]

        #Assign external nodes/tips to OUi
        trdata$dat$regimes[external.nodes]<-paste("OU",sep="",i+1)
        #Assign internal nodes to OUi
        internal.nodes.regimes[internal.nodes-n.tips]<-paste("OU",sep="",i+1)
        #Assign node where shift occurs to OUi
        internal.nodes.regimes[rep.regimes-n.tips]<-paste("OU",sep="",i+1)
        }
      }
    }

  trdata$phy$node.label<-internal.nodes.regimes
  reg.colors<-ggsci::pal_npg(palette=c("nrc"),alpha=1)(length(regimes)+1)
  #reg.colors<-microViz::distinct_palette(n=length(regimes)+1)

  #Combine external coding and internal coding to plot tree with colored shifts
  regimes.total<-c(trdata$dat$regimes,internal.nodes.regimes)
  edge.regimes <- factor(regimes.total[trdata$phy$edge[,2]])
  print(edge.regimes)
  print(reg.colors)
  plot(trdata$phy,edge.color = reg.colors[edge.regimes], edge.width = 1, cex = 0.2)
  #return(list(trdata,internal.nodes.regimes))
  return(trdata)
}
