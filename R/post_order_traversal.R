library(ape)

# Create a simple 5-tip tree
set.seed(123) # For reproducibility
n_tips<-5
phy <- rtree(n_tips)

# Plot the tree to visualize it (optional, but helpful)
plot(phy, show.node.label = TRUE) # show.node.label helps see node IDs
nodelabels(frame = "n", cex = 0.8, col = "blue") # Add node numbers
tiplabels(frame = "n", cex = 0.8, col = "red") # Add tip numbers
num_nodes<-phy$Nnode+5
phy$edge

parent_of_node <- rep(0, num_nodes) #Asign all 0s first, equal to number of nodes
#For each child ID, will store the parent of that child, parent_of_node[child_id] = parent_id
children_of_node_list <- vector("list", num_nodes) #Create vector of list equal to the number of nodes - this will store the children of each parent node - given a parent node id, will return all child node ids
#children_of_node_list[[parent_id]] = c(child1_id, child2_id, ...)

for(i in 1:nrow(phy$edge)){ #For all branches - parent to child nodes in order of node numbers -
  parent <- phy$edge[i, 1] #Get parent node
  child <-phy$edge[i,2] #Get child node
  parent_of_node[child] <- parent #For vector of parent ids, for element = number of child, place parent id in it - returns vector with parent ID for each child ID
  children_of_node_list[[parent]] <- c(children_of_node_list[[parent]], child) #For list of child ids, for the parent node, add on all children that come from that parent node - returns list with children IDs for each parent ID - NULL for tips. Only returns immediate children of parent
}


#Left - right - root
stack1 <-root_node_id
stack2 <- integer(0)

while(length(stack1)>0){ # LIFO (Last-In, First-Out) stack; https://www.naukri.com/code360/library/postorder-traversal-of-binary-tree
  current<-tail(stack1,n=1)
  stack1<-stack1[-length(stack1)]
  stack2<-c(current,stack2)
  stack1<-c(stack1,children_of_node_list[[current]])
}

