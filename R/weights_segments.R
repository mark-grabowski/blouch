#' weights_segments - for internal Blouch use
#' For individual lineage, determine the weighting of each segment
#' @param a Rate parameter for OU model
#' @param lineage Individual lineage regime values
#'
#' @return For individual lineage, determine the weighting of each segment
#' @export
#'
weights_segments <- function(a, lineage){#For individual lineage, determine the weighting of each segment
  #t_beginning and t_end are both vectors, and subtracting them from each other lines up the beginning and end of one segment
  #because if tge tail/head procedure above
  res <- c(exp(-a * lineage$t_beginning) - exp(-a * lineage$t_end),
           exp(-a * lineage$times[1]))
  return(res)
}
