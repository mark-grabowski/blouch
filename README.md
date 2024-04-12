
<!-- README.md is generated from README.Rmd. Please edit that file -->

# *Blouch*

<!-- badges: start -->
<!-- badges: end -->

*Blouch*: Bayesian Linear Ornstein-Uhlenbeck models for Comparative
Hypotheses fits allometric and adaptive models of continuous trait
evolution in a Bayesian framework based on fixed or continuous
predictors and incorporates measurement error. In addition to assigning
biologically meaningful priors when compared to non-Bayesian approaches,
*Blouch* includes new implementations of Ornstein-Ulenbeck models
including allowing for varying effects (varying intercepts and varying
slopes), multilevel modeling, and non-centered models.

While the front-end component of *Blouch* is written in R (R Core Team,
2023), the nuts and bolts are written in the language Stan (Carpenter et
al., 2017), which allows estimation of Bayesian models using Markov
chain Monte Carlo (MCMC) methods based on the Hamilton Monte Carlo
sampler.

## Getting Started

If you are just getting started with *Blouch* I recommend starting with
the Empirical and Simulation Example articles available on the package website. The
other articles are abbreviated versions of this example showing the
various models implemented by *Blouch* - most of the steps of a
preliminary analysis will be the same.

*Blouch* is based on an article currently in review:

- Grabowski, M (in revision). *Blouch*: Bayesian Linear
  Ornstein-Uhlenbeck models for Comparative Hypotheses

## Instalation Instructions

To install the R and Stan functions associated with Blouch from github,
first install the package devtools:

``` r
install.packages("devtools", repos = "https://cran.ma.imperial.ac.uk/")
library(devtools)
```

Then install Blouch

``` r
devtools::install_github("mark-grabowski/blouch")
library(blouch)
```

## Documentation

For a complete walkthrough of the package visit the package website:
<a href="https://mark-grabowski.github.io/blouch/" title="here.">here</a>.
Package development is ongoing at the blouch-project github page:
<a href="https://mark-grabowski.github.io/blouch-project/" title="here.">here</a>.

## References

Carpenter, B., A. Gelman, M. D. Hoffman, D. Lee, B. Goodrich, M.
Betancourt, M. Brubaker, J. Guo, P. Li, and A. Riddell. 2017. Stan: A
Probabilistic Programming Language. Journal of Statistical Software
76:1–32.

R Core Team. 2023. R: A language and environment for statistical
computing.
