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
#' @param obs_gamma observed gamma value to fit on
#' @param obs_colless observed colless value to fit on
#' @param obs_num_lin observed number of lineages to fit on
#' @param limiting_accept_rate limiting accept rate
#' @param sim_number sim number for Thijs
#' @param bd_lambda lambda estimate using DDD::bd_ML
#' @param bd_mu mu estimate using DDD::bd_ML
#' @return tibble
#' @export
#' @rawNamespace import(Rcpp)
#' @rawNamespace useDynLib(pbdABC, .registration = TRUE)
perform_abc_serial <- function(num_particles,
                            num_iterations,
                            crown_age,
                            min_lin,
                            max_lin,
                            lower = c(-3, -3, -3, -3, -6),
                            upper = c(3, 3, 3, 3, 6),
                            means = c(-1, -1, -1, -1, -1),
                            obs_gamma,
                            obs_colless,
                            obs_num_lin,
                            limiting_accept_rate = 1e-8,
                            sim_number,
                            bd_lambda = 0,
                            bd_mu = 0) {

  res <- perform_abc_rcpp(num_particles,
                          num_iterations,
                          crown_age,
                          min_lin,
                          max_lin,
                          lower,
                          upper,
                          means,
                          obs_gamma,
                          obs_colless,
                          obs_num_lin)

  #colnames(res) <- c("iter", "lambda0", "lambda1", "mu0", "mu1", "compl_rate",
  #                   "gamma", "colless", "num_lin", "weight")
  #res <- tibble::as_tibble(res)

  return(res)
}
