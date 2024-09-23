################################################################################
#
# @brief Estimate the likelihood of a given tree, provided a likelihood
#'        function, using a Monte Carlo Markov Chain
#'        Lifted from the nltt package
#'
# @date Last modified: 2014-20-09
# @author Thijs Janzen
# @since 2014-20-09, version 1.0
#'
#' @param    phy                   phylo       Vector of weights
#' @param    likelihood_function   function    Function that calculates the
#'                                             likelihood of our diversification
#'                                             model, given the tree.
#'                                             function should be of the format
#'                                             function(parameters, phy).
#' @param    parameters            vector      Initial parameters to start
#'                                             the chain.
#' @param    logtransforms         scalar      Whether to perform jumps on
#'                                             logtransformed parameters (TRUE)
#'                                             or not (FALSE)
#' @param    iterations            scalar      Length of the chain
#' @param    burnin                scalar      Length of the burnin, default is
#'                                             30\% of iterations
#' @param    thinning              scalar      Size of thinning, default = 1
#' @param    sigma                 scalar      Standard deviation of the jumping
#'                                             distribution, which is
#'                                             N(0, sigma).
#' @param    sim_number           sim number
#' @return                         mcmc        An MCMC object, as used by the
#'                                             package "coda".
#'
################################################################################
#' @export
mcmc_pbd <- function(phy, likelihood_function,
                     parameters, logtransforms, iterations,
                     burnin = round(iterations / 3),
                     thinning = 1,
                     sigma = 1,
                     sim_number
) {

  #check data type of phy
  if (!inherits(phy, "phylo")) {
    # Just checking
    stop("mcmc: ",
         "phy must be of class 'phylo', ",
         "but was of type '", class(phy), "' instead")
  }

  # create a list for the samples & reserve memory for the chain
  chain <- array(dim = c(floor(iterations / thinning) + 1,
                         length(parameters)))

  for (j in seq_along(parameters)) {
    if (parameters[j] < 0) {
      #Just checking
      stop("mcmc: ",
           "initial parameter values have to be above zero\n",
           "but one was ", parameters[j], " instead")
    }
  }
  # pre-compute current posterior probability
  pp <- likelihood_function(parameters, phy)

  cat("\nGenerating Chain\n")
  cat("0--------25--------50--------75--------100\n")
  cat("*")
  utils::flush.console()
  print_frequency <- 20

  file_name <- paste0("mcmc_", sim_number, ".txt")

  for (i in seq_len(burnin + iterations)) {
    #propose new values
    for (j in seq_along(parameters)) {
      if (logtransforms[j] == TRUE) {
        if (parameters[j] == 0) {
          stop("Cannot propose new value for a parameter with value 0.0.")
        }

        eta           <- log(parameters[j])
        new_eta       <- eta + stats::rnorm(1, 0, sigma)
        new_val       <- exp(new_eta)
        # calculate the Hastings ratio
        hr            <- log(new_val / parameters[j])
        parameters[j] <- new_val
        new_pp        <- likelihood_function(parameters, phy)

        #accept or reject
        if (is.finite(new_pp) &&
            is.finite(hr) &&
            new_pp - pp + hr > log(stats::runif(1, 0, 1))) {
          pp <- new_pp
        } else {
          parameters[j] <- exp(eta)
        }
      } else {

        eta           <- parameters[j]
        new_val       <- eta + stats::rnorm(1, 0, sigma)
        #calculate the Hastings ratio
        hr            <- 0.0
        parameters[j] <- new_val

        if (parameters[j] >= 0 & parameters[1] > 0) {
          new_pp        <- likelihood_function(parameters, phy)

          #accept or reject
          if (is.finite(new_pp) &&
              is.finite(hr) &&
              new_pp - pp + hr > log(stats::runif(1, 0, 1))) {
            pp <- new_pp
          } else {
            parameters[j] <- eta
          }
        } else {
          parameters[j] <- eta
        }
      }
    }

    # sample the parameter
    if (i >= burnin) {
      if ((i) %% ((iterations - burnin) / print_frequency) == 0) {
        cat("**")
        utils::flush.console()
      }
      if ((i - burnin) %% thinning == 0) {
        chain[(i - burnin) / thinning + 1, ] <- parameters
        cat(parameters, "\n", file = file_name, append = TRUE)
      }
    }
  }
  cat("\nFinished MCMC.\n")
  #return a mcmc object, used by coda to plot
  return(coda::as.mcmc(chain))
}
