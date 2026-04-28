#Four regimes with one adaptive trait and multiple slopes per optima but single alpha parameter
set.seed(10) #Set sequence of random numbers for reproducability
N<-3 #Number of species

phy <- ape::keep.tip(tree.10K,sample(tree.10K$tip.label)[1:N])
phy<-ape::multi2di(phy) #Collapse or resolve multichotomies in phylogenetic trees.

l.tree<-max(ape::branching.times(phy)) ## Rescale tree to height 1
phy$edge.length<-phy$edge.length/l.tree
plot(phy)
ape::nodelabels(frame="none",adj=c(1.1,-0.4))
ape::tiplabels()



calc_path_matrix<-function(phy){
  N.tips<-length(phy$tip.label)
  N.branches <- nrow(phy$edge)
  path_mat <- matrix(0, nrow = N.tips, ncol = N.branches) # Set up path matrix
  for (i in 1:(N.tips)){
    path_to_tip <- ape::nodepath(phy, from = N.tips + 1, to = i) # For each tip, determine its path from root to tip
    for(j in 1:(length(path_to_tip)-1)){ # For each path, determine which branches fall on that path, -1 for walking through path
      parent <- path_to_tip[j] # Parent of path
      child <- path_to_tip[j+1] # Child of path
      edge_idx <- which(phy$edge[,1] == parent & phy$edge[,2] == child) # Which edge matches that parent-child pair?
      path_mat[i, edge_idx ] <- 1 # Assign 1's in path matrix [N.tips X N.braches/edges] that have that paid
    }
  }
  return(path_mat)
}

calc_branch_matrix<-function(phy){
  N.tips<-length(phy$tip.label)
  N.branches <- nrow(phy$edge)
  branch_mat <- matrix(0, nrow = N.tips, ncol = N.branches) # Set up path matrix
  for (i in 1:(N.tips)){
    path_to_tip <- ape::nodepath(phy, from = N.tips + 1, to = i) # For each tip, determine its path from root to tip
    for(j in 1:(length(path_to_tip)-1)){ # For each path, determine which branches fall on that path, -1 for walking through path
      parent <- path_to_tip[j] # Parent of path
      child <- path_to_tip[j+1] # Child of path
      edge_idx <- which(phy$edge[,1] == parent & phy$edge[,2] == child) # Which edge matches that parent-child pair?
      branch_mat[i, edge_idx ] <- phy$edge.length[edge_idx] # Assign 1's in path matrix [N.tips X N.braches/edges] that have that paid
    }
  }
  return(branch_mat)
}


path_mat <- calc_path_matrix(phy)
#branch_mat <- calc_branch_matrix(phy)
branch_vec <- phy$edge.length
times.MRCA<-stats::cophenetic(phy)



time_to_end_node <- ape::node.depth.edgelength(phy)
child.nodes <- phy$edge

# End of each lineage is tip, want sum of path from root to tip


N.tips<-length(phy$tip.label)
N.branches <- nrow(phy$edge)

time_to_tip <- rep.int(0, N.tips)

for (i in 1:N){
  for (j in 1:M){
    if (path_mat[i, j] ==1 ){
      weighted_path[i, j] =

    }

  }

}

for(i in 1:N.branches){ #
  run.sum <- 0
  for(j in 1:N.branches){
    if(path_mat[i, j] == 1)
      {
      run.sum <- run.sum + phy$edge.length[j]
      }
    time_to_tip[i] <- run.sum
}
}

calc_time_to_tips_matrix <- function(phy, path_mat)
  N.tips<-length(phy$tip.label)
  N.branches <- nrow(phy$edge)
  time_mat <- matrix(0, nrow = N.tips, ncol = N.branches) # Set up time matrix
  children <- phy$edge[,2]
  branch_lengths <- rep.int(0, N.branches) # Default = 0 branch lengths
  for (i in 1:(N.tips)){
    for(j in 1:(N.branches)){ # For each path, determine which branches fall on that path, -1 for walking through path
      if(path_mat[i, j] == 1 ){

      }

      child <- path_to_tip[j+1] # Child of path
      edge_idx <- which(phy$edge[,1] == parent & phy$edge[,2] == child) # Which edge matches that parent-child pair?
      path_mat[i, edge_idx ] <- 1 # Assign 1's in path matrix [N.tips X N.braches/edges] that have that paid
    }
  }
  return(path_mat)
}





}
for (i in 1:N.tips){



}


blouch.reg.shifts.prep<-function(phy,Y,Y_error,hl.prior,vy.prior,optima.prior){
    #Data should be send in treeplyr format
    phy<-trdata$phy
    dat<-data.frame(trdata$dat)
    N.tips<-length(trdata$phy$tip.label)

    ############################################################################################################
    Dmat<-stats::cophenetic(phy) # Time separating tips
    time_to_end_node <- ape::node.depth.edgelength(phy)
    path_mat <- calc_path_matrix(phy) #Calculate path matrix
    #
    ############################################################################################################
    #print(as.vector(t(dat[Y_error])))
    dat<-list(N=N,
              Y_obs=as.vector(t(dat[Y])),
              Y_error=as.vector(t(dat[Y_error])),
              path_mat,
              Dmat=Dmat,
              hl_prior=as.vector(hl.prior),
              vy_prior=vy.prior,
              optima_prior=as.vector(optima.prior))

    return(dat)
}

blouch.reg.shifts.prep()


