#' function to do abc in parallel
#' @param num_particles number of particles
#' @param num_iterations number of iterations
#' @param crown_age crown age
#' @param min_lin minimum number of lineages
#' @param max_lin maximum number of lineages
#' @param lambdas vector of lambdas for exponential priors (5),
#' order: l0, mu0, l1, mu1, compl_rate
#' @param obs_gamma observed gamma value to fit on
#' @param obs_colless observed colless value to fit on
#' @param obs_num_lin observed number of lineages to fit on
#' @param num_threads number of threads
#' @return matrix
#' @export
#' @rawNamespace import(Rcpp)
#' @rawNamespace useDynLib(pbdABC, .registration = TRUE)
test_abc_par <- function(num_particles,
                            num_iterations,
                            crown_age,
                            min_lin,
                            max_lin,
                            lambdas,
                            obs_gamma,
                            obs_colless,
                            obs_num_lin,
                            num_threads = 1) {
  RcppParallel::setThreadOptions(numThreads = num_threads)
  res <- test_abc_rcpp_par(num_particles,
                              num_iterations,
                              crown_age,
                              min_lin,
                              max_lin,
                              lambdas,
                              obs_gamma,
                              obs_colless,
                              obs_num_lin)
  return(res)
}
