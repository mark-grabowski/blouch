#' parent function - Returns parent node of offspring node given node number
#'
#' @param phy phylogeny in phytools format
#' @param x node number
#'
#' @return value Returns parent node of offspring node given node number
#' @export
#'
parent <- function(phy, x){ #Returns parent node of offspring node given node number
  m <- which(phy$edge[, 2] == x)
  return(phy$edge[m, 1])
}
