#' function to do abc in parallel
#' @param num_particles number of particles
#' @param num_iterations number of iterations
#' @param crown_age crown age
#' @param min_lin minimum number of lineages
#' @param max_lin maximum number of lineages
#' @param lower lower prior limits
#' @param upper upper prior limits
#' @param obs_gamma observed gamma value to fit on
#' @param obs_colless observed colless value to fit on
#' @param obs_num_lin observed number of lineages to fit on
#' @param num_threads number of threads
#' @param sim_number sim number for Thijs
#' @return matrix
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
                      sim_number) {
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
                              sim_number)
  return(res)
}
