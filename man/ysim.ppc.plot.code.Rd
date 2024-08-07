% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create.plots.R
\name{ysim.ppc.plot.code}
\alias{ysim.ppc.plot.code}
\title{ysim.ppc.plot.code - Create plot for Prior and Poserior Predictive Checks}
\usage{
ysim.ppc.plot.code(dat, post, row.nums)
}
\arguments{
\item{dat}{Data formatted for Blouch models by blouch.prep functions}

\item{post}{Posterior distribution of stanfit class}

\item{row.nums}{Rows to be sampled from distribution}
}
\value{
Plots in ggplot2 format
}
\description{
ysim.ppc.plot.code - Create plot for Prior and Poserior Predictive Checks
}
