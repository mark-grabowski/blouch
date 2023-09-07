#' weight.matrix - for internal Blouch use
#' Wrapper to apply weights_regimes to each lineage
#' @param phy phylogeny in phytools format
#' @param a OU rate parameter
#' @param lineages Vector of regime values
#'
#' @return weights for each lineage
#' @export
#'
weight.matrix <- function(phy, a, lineages){ #Wrapper to apply weights_regimes to each lineage
  res <- t(vapply(lineages, function(x) weights_regimes(a, x),
                  FUN.VALUE = numeric(length(lineages[[1]]$which.regimes))) ## Specify type of output
  )

  rownames(res) <- phy$tip.label
  return(res)
}
