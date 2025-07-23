root_trace_redux<-function(trdata,N){ # Needs a formatted trdata file from set.converge.regimes

  # Extract necessary components from trdata
  phy<-trdata$phy # Store phylogeny
  tip_regimes<-trdata$dat$regimes # Store regimes at tips (this is the regime of the branch leading *to* the tip)

  # Prepare edge table for efficient lookup
  edge_table<-phy$edge # Parent, child tip/node
  edge_table<-cbind(edge_table,phy$edge.length) # Bind numbers and lengths of edges together
  edge_table<-data.frame(edge_table)
  names(edge_table)<-c("parent","child","edge_length")

  root_num<-N+1 # Root node is N+1 in ape convention for an N-tip tree

  # Initialize lists to store results for each tip's path
  store_path_list<-list()       # Sequence of node numbers from root to tip
  store_reg_list<-list()        # Sequence of regimes for each branch segment on the path
  store_edge_length_list<-list()# Sequence of branch lengths for each segment on the path

  # Extract internal node regimes (as assigned by set.converge.regimes)
  # Based on set.converge.regimes, node_regimes[k-N] is the regime of the branch *descending from* internal node k.
  node_regimes<-phy$node.label

  # Loop through each tip to trace its path back to the root
  for(i in 1:N){
    current_child <- i # Start tracing from the current tip
    path_nodes <- c(current_child) # Initialize path with the tip node itself

    # Initialize lists to build up regimes and lengths for this specific path
    # We'll prepend elements as we trace back up the tree, then these lists will be in root-to-tip order.
    branch_regimes <- c()
    branch_lengths <- c()

    # Loop to trace back from the current child (initially a tip) to the root
    while(TRUE){
      # Find the edge that leads *to* the current_child
      edge_idx <- which(edge_table$child == current_child)

      # Break if current_child is the root (it won't have a parent edge in edge_table)
      # This effectively handles the root as the end of the path.
      if (length(edge_idx) == 0) {
        break
      }

      current_parent <- edge_table$parent[edge_idx]
      current_edge_length <- edge_table$edge_length[edge_idx]

      # Prepend the branch length to maintain root-to-tip order
      branch_lengths <- c(current_edge_length, branch_lengths)

      # Determine the regime of the branch (current_parent -> current_child)
      # If current_child is a tip, its regime comes from trdata$dat$regimes (tip_regimes)
      # If current_child is an internal node, its regime is associated with its parent
      # Based on set.converge.regimes, node_regimes[parent - N] stores the regime of the branch *descending* from that parent.
      # So, the regime of the branch (current_parent -> current_child) is stored at current_parent's node label.
      # (This mapping needs to be accurate: internal node labels refer to branches *descending* from them).

      if (current_child <= N) { # If the current_child is a tip
        regime_to_prepend <- tip_regimes[current_child]
      } else { # If the current_child is an internal node
        # The regime of the branch (current_parent -> current_child) is the regime assigned to current_parent's label.
        # node_regimes is indexed from 1 for node N+1, so parent - N corrects the index.
        regime_to_prepend <- node_regimes[current_parent - N]
      }
      branch_regimes <- c(regime_to_prepend, branch_regimes)


      # Prepend the parent node to the path
      path_nodes <- c(current_parent, path_nodes)

      # Move up the tree: the current parent becomes the new child for the next iteration
      current_child <- current_parent

      # If we just added the root, we're done with this path
      if (current_child == root_num) {
        break
      }
    }

    # Store the complete path, regimes, and branch lengths for the current tip
    # These lists are now in root-to-tip order.
    store_path_list[[length(store_path_list)+1]]<-path_nodes
    store_reg_list[[length(store_reg_list)+1]]<-branch_regimes
    store_edge_length_list[[length(store_edge_length_list)+1]]<-branch_lengths
  }

  # Return all collected paths
  return(list(path_nodes = store_path_list,
              regimes = store_reg_list,
              branch_lengths = store_edge_length_list))
}
