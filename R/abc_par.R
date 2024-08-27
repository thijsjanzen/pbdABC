#' function to do abc in parallel
#' @param num_particles number of particles
#' @param num_iterations number of iterations
#' @param crown_age crown age
#' @param min_lin minimum number of lineages
#' @param max_lin maximum number of lineages
#' @param lower lower prior limits, ordering: l0, l1, m0, m1, compl_rate
#' @param upper upper prior limits, ordering: l0, l1, m0, m1, compl_rate
#' @param obs_gamma observed gamma value to fit on
#' @param obs_colless observed colless value to fit on
#' @param obs_num_lin observed number of lineages to fit on
#' @param num_threads number of threads
#' @param limiting_accept_rate limiting accept rate
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
                      lower,
                      upper,
                      obs_gamma,
                      obs_colless,
                      obs_num_lin,
                      num_threads = 1,
                      limiting_accept_rate = 1e-8,
                      sim_number,
                      bd_lambda = 0,
                      bd_mu = 0) {
  RcppParallel::setThreadOptions(numThreads = num_threads)
  res <- perform_abc_rcpp_par(num_particles,
                              num_iterations,
                              crown_age,
                              min_lin,
                              max_lin,
                              lower,
                              upper,
                              obs_gamma,
                              obs_colless,
                              obs_num_lin,
                              sim_number,
                              bd_lambda,
                              bd_mu,
                              limiting_accept_rate)

  colnames(res) <- c("iter", "lambda0", "lambda1", "mu0", "mu1", "compl_rate",
                     "gamma", "colless", "num_lin")
  res <- tibble::as_tibble(res)

  return(res)
}
