#' The 'blouch' package.
#'
#' @description Bayesian Linear Ornstein-Uhlenbeck models for Comparative Hypotheses (BLOUCH) fits adaptive models of continuous trait evolution in a Bayesian framework based on categorical or continuous predictors, and incorporates measurement error.
#'
#' @docType package
#' @name blouch-package
#' @aliases blouch
#' @import methods
#' @import Rcpp
#' @import parallel
#' @import StanHeaders
#' @importFrom geiger ratematrix
#' @importFrom utils head tail
#' @importFrom rstan sampling
#' @importFrom RcppParallel RcppParallelLibs
#' @useDynLib blouch, .registration = TRUE
#'
#' @references
#' Stan Development Team (2023). RStan: the R interface to Stan. R package version 2.26.22. https://mc-stan.org
#'
NULL
