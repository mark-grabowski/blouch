set.converge.regimes<-function(trdata,regimes){
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

  trdata$dat$regimes<-"OU1"
  internal.nodes.regimes<-rep("OU1",n.internal.nodes)

  for(i in 1:(length(regimes))){
    # print(i) # Consider removing print statements for performance
    # Corrected loop structure: The second inner 'j' loop was redundant and should be merged
    for(j in 1:length(regimes[[i]])){
      # This 'if' block handles shifts directly on tip branches (ancestor is an internal node, descendant is a tip)
      if(regimes[[i]][[j]] <= n.tips){ # shift on single branch leading to a tip
        trdata$dat$regimes[regimes[[i]][[j]]] <- paste("OU", sep = "", i + 1)
      }

      rep.regimes <- regimes[[i]][[j]] # The node where the shift occurs (can be tip or internal node)
      saved.decendants <- getDescendants(trdata$phy, node = rep.regimes)
      external.nodes <- saved.decendants[saved.decendants <= n.tips]
      internal.nodes <- saved.decendants[saved.decendants > n.tips]

      # Assign external nodes/tips to OUi (branches leading to these tips)
      trdata$dat$regimes[external.nodes] <- paste("OU", sep = "", i + 1)
      # Assign internal nodes to OUi (branches descending from these internal nodes)
      internal.nodes.regimes[internal.nodes - n.tips] <- paste("OU", sep = "", i + 1)
      # Assign the node where shift occurs to OUi (the branch descending from this shift node)
      internal.nodes.regimes[rep.regimes - n.tips] <- paste("OU", sep = "", i + 1)
    }
  }

  trdata$phy$node.label<-internal.nodes.regimes
  # The plotting part is fine for visualization
  reg.colors<-ggsci::pal_npg(palette=c("nrc"),alpha=1)(length(regimes)+1)
  regimes.total<-c(trdata$dat$regimes,internal.nodes.regimes)
  edge.regimes <- factor(regimes.total[trdata$phy$edge[,2]]) # Assigns regime based on the CHILD node's regime
  print(edge.regimes) # Consider removing print statements for performance
  print(reg.colors)   # Consider removing print statements for performance
  plot(trdata$phy,edge.color = reg.colors[edge.regimes], edge.width = 1, cex = 0.2)
  return(trdata)
}
