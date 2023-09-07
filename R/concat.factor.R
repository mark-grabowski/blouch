#' concat.factor - for internal Blouch use
#' Thanks to user "snaut" at stackoverflow, http://stackoverflow.com/users/1999873/snaut
#' @param ... vector of factors
#'
#' @return factor
#' @export
#'
concat.factor <- function(...){
  as.factor(do.call(c, lapply(list(...), as.character)))
}
