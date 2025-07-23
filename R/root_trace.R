root_trace<-function(trdata,N){ #Needs a formatted trdata file from set.converge.regimes
  phy<-trdata$phy #Store phylogeny
  tip_regimes<-trdata$dat$regimes #Store regimes at tips
  edge_table<-phy$edge #Parent, child tip/node
  edge_table<-cbind(edge_table,phy$edge.length) #Bind numbers and lengths of edges together into a table
  edge_table<-data.frame(edge_table)
  names(edge_table)<-c("parent","child","edge_length")
  root_num<-N+1 #Root is N+1 always
  store_path_list<-list()
  store_reg_list<-list()
  store_edge_length_list<-list()
  node_regimes<-phy$node.label #Regimes of the nodes in the order of the node numbers
  for(i in 1:N){
    edge_length<-NULL
    edge_num<-match(i,edge_table$child) #Match the tip number with the row number of the edge table
    child <- edge_table$child[edge_num] #Extract child and parent node number
    parent <- edge_table$parent[edge_num]
    store_path<-c(parent,child) #Store in store_path variable - both
    store_reg<-c(node_regimes[parent-N],tip_regimes[i]) #Store parent and tip regime data
    store_edge_length<-edge_table$edge_length[edge_num]
    while(parent != root_num){ #Run until the parent number is the root number
      edge_num<-match(parent,edge_table$child) #Find the edge number of the old parent that is the new child
      child <- edge_table$child[edge_num] #Store parent and child node number
      parent <- edge_table$parent[edge_num]
      store_path<-c(parent,store_path) #Store parent node number in store_path variable
      store_reg<-c(node_regimes[parent-N],store_reg) #Store parent regime in store_reg
      store_edge_length<-c(edge_table$edge_length[edge_num],store_edge_length) #Store edge_length
    }
    store_path_list[[length(store_path_list)+1]]<-store_path #Starts with root and works it way to the end
    store_reg_list[[length(store_reg_list)+1]]<-store_reg #Starts with root and works it way to the end
    store_edge_length_list[[length(store_edge_length_list)+1]]<-store_edge_length #Starts with root and works it way to the end

  }

return(list(store_path_list,store_reg_list,store_edge_length_list))
}
