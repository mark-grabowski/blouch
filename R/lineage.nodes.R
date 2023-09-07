#' lineage.nodes function for internal Blouch use
#' Given a certain node, return the list of all parent nodes back to the root of the tree
#' @param phy phylogeny in phytools format
#' @param x node of interest
#'
#' @return list Given a certain node, return the list of all parent nodes back to the root of the tree
#' @export
#'
lineage.nodes <- function(phy, x){ #Given a certain node, return the list of all parent nodes back to the root of the tree
  k <- x #Save x in k
  N <- length(phy$tip.label) #total number of tips on tree
  while(x != N + 1){ #while node is not equal to number of tips +1  - starting node - 51 here
    k <- c(k, parent(phy, x)) #Return node at beginning of edge
    x <- tail(k, n = 1) #x is assigned value at end of k, so end of the list of beginning nodes
    #50->99->89->51 0- tracing lineage back by nodes
  }
  return(k)
}
