#' function to do abc in parallel
#' @param num_particles number of particles
#' @param num_iterations number of iterations
#' @param crown_age crown age
#' @param min_lin minimum number of lineages
#' @param max_lin maximum number of lineages
#' @param lower lower prior limits, ordering: l0, l1, m0, m1, compl_rate.
#' Given are the log10 values of the prior, e.g. -5 = 10^-5.
#' @param upper upper prior limits, ordering: l0, l1, m0, m1, compl_rate.
#' Given are the log10 values of the prior, e.g. 5 = 10^5.
#' @param means mean prior values, when using an exponential prior. Values
#' should be > 0. Default are negative values, in which case the log-uniform
#' prior is used.
#' @param use_inv_prior use of an inverse prior on completion rate?
#' @param obs_gamma observed gamma value to fit on
#' @param obs_colless observed colless value to fit on
#' @param obs_num_lin observed number of lineages to fit on
#' @param num_threads number of threads
#' @param limiting_accept_rate stop when the acceptance rate of a generation
#' drops below this value
#' @param threshold_quantile NA (default) uses the fixed threshold schedule
#' 10 * exp(-0.5 * (i - 1)). A value between 0 and 1 uses an adaptive
#' threshold: each generation accepts particles whose distance is at most
#' this quantile of the distances of the previous generation (0.5 is a
#' common choice).
#' @param sim_number sim number for Thijs
#' @param bd_lambda lambda estimate using DDD::bd_ML
#' @param bd_mu mu estimate using DDD::bd_ML
#' @return tibble
#' @export
#' @rawNamespace import(Rcpp)
#' @rawNamespace useDynLib(pbdABC, .registration = TRUE)
perform_abc_par <- function(num_particles,
                                 num_iterations,
                        crown_age,
                      min_lin,
                      max_lin,
                      lower = c(-3, -3, -3, -3, -6),
                      upper = c(3, 3, 3, 3, 6),
                      means = c(-1, -1, -1, -1, -1),
                      use_inv_prior = FALSE,
                      obs_gamma,
                      obs_colless,
                      obs_num_lin,
                      num_threads = 1,
                      limiting_accept_rate = 1e-12,
                      sim_number,
                      bd_lambda = 0,
                      bd_mu = 0,
                      threshold_quantile = NA) {
  if (!is.na(threshold_quantile) &&
      (threshold_quantile <= 0 || threshold_quantile >= 1)) {
    stop("threshold_quantile should be NA or between 0 and 1")
  }
  RcppParallel::setThreadOptions(numThreads = num_threads)
  res <- perform_abc_rcpp_par(num_particles,
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
                              obs_num_lin,
                              sim_number,
                              bd_lambda,
                              bd_mu,
                              limiting_accept_rate,
                              ifelse(is.na(threshold_quantile), -1,
                                     threshold_quantile))

  colnames(res) <- c("iter", "lambda0", "lambda1", "mu0", "mu1", "compl_rate",
                     "gamma", "colless", "num_lin", "weight", "dist")
  res <- tibble::as_tibble(res)

  return(res)
}
