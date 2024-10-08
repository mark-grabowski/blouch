% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/set.converge.regimes.R
\name{set.converge.regimes}
\alias{set.converge.regimes}
\title{set.converge.regimes - function to assign regimes on a phylogeny}
\usage{
set.converge.regimes(trdata, regimes)
}
\arguments{
\item{trdata}{treeplyr format file}

\item{regimes}{node numbers for where the regimes shift}
}
\value{
treeplyr format file with new column in dat for tip regimes and internal regime assignments in node.labels
}
\description{
set.converge.regimes - function to assign regimes on a phylogeny
}
