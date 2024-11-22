% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{test_abc_rcpp_par}
\alias{test_abc_rcpp_par}
\title{function to do abc using rcpp}
\usage{
test_abc_rcpp_par(
  num_particles,
  num_iterations,
  crown_age,
  min_lin,
  max_lin,
  lower,
  upper,
  means,
  use_inv_prior,
  obs_gamma,
  obs_colless,
  obs_num_lin
)
}
\arguments{
\item{num_particles}{number of particles}

\item{num_iterations}{number of iterations}

\item{crown_age}{crown age}

\item{min_lin}{minimum number of lineages from the prior}

\item{max_lin}{maximum number oflineages from the prior}

\item{lower}{minimum values of log-uniform prior (e.g. -3 = 10^-3), ordering:
l0, l1, mu0, mu1, compl_rate}

\item{upper}{upper values of log-uniform prior (e.g. 3 = 10^3)//' @param obs_gamma observed gamma value to fit on}

\item{use_inv_prior}{use inverse prior}

\item{obs_gamma}{observed gamma value to fit on}

\item{obs_colless}{observed colless value to fit on}

\item{obs_num_lin}{observed number of lineages to fit on}
}
\description{
Fit to the data is assessed by the sum of differences for
[gamma, colless, num_lineages], rescaled by the observed value, e.g.
((O-E)^2)/E, where O is the value of the proposed simulation and E is the
value of the empirical data (e.g. obs_gamma, obs_colless or obs_num_lineages).
The acceptance threshold diminishes exponentially.
}
