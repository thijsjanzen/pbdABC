# Plotting functions for the output of perform_abc_par().

abc_param_names <- c("lambda0", "lambda1", "mu0", "mu1", "compl_rate")

abc_param_labels <- c(lambda0    = "lambda0: initiation (good)",
                      lambda1    = "lambda1: initiation (incipient)",
                      mu0        = "mu0: extinction (good)",
                      mu1        = "mu1: extinction (incipient)",
                      compl_rate = "compl_rate: completion")

utils::globalVariables(c("value", "weight", "parameter", "iter",
                         "density", "ndensity", "x", "generation"))

# long format: one row per particle and parameter, weights normalised
# within each generation
abc_to_long <- function(res, iterations = NULL) {
  res <- as.data.frame(res)
  needed <- c("iter", abc_param_names, "weight")
  missing <- setdiff(needed, colnames(res))
  if (length(missing) > 0) {
    stop("res should be the output of perform_abc_par; missing columns: ",
         paste(missing, collapse = ", "))
  }
  if (!is.null(iterations)) res <- res[res$iter %in% iterations, ]
  if (nrow(res) == 0) stop("no particles to plot")

  res$weight <- res$weight / stats::ave(res$weight, res$iter, FUN = sum)

  long <- do.call(rbind, lapply(abc_param_names, function(p) {
    data.frame(iter = res$iter, parameter = p, value = res[[p]],
               weight = res$weight)
  }))
  long$parameter <- factor(abc_param_labels[long$parameter],
                           levels = abc_param_labels)
  long
}

# Weighted kernel density of each parameter in each generation, on the
# log10 scale. The bandwidth follows Silverman's rule with the weighted
# standard deviation and the effective sample size (stats::density's own
# bandwidth selectors ignore the weights). If log10 prior bounds are given,
# the density is restricted to them, using reflection at the bounds so that
# no probability mass is lost near the edges.
abc_weighted_density <- function(long, lower = NULL, upper = NULL, n = 512) {
  groups <- split(long, list(long$iter, long$parameter), drop = TRUE)
  out <- lapply(groups, function(g) {
    x <- log10(g$value)
    w <- g$weight / sum(g$weight)
    ess <- 1 / sum(w^2)
    m <- sum(w * x)
    s <- sqrt(sum(w * (x - m)^2))
    bw <- 0.9 * s * ess^(-1 / 5)
    if (!is.finite(bw) || bw <= 0) bw <- 0.01  # all particles identical
    if (is.null(lower) || is.null(upper)) {
      d <- stats::density(x, weights = w, bw = bw, n = n)
      y <- d$y
    } else {
      j <- as.integer(g$parameter[1])  # factor levels follow abc_param_names
      lo <- lower[j]
      up <- upper[j]
      d <- stats::density(c(x, 2 * lo - x, 2 * up - x),
                          weights = rep(w, 3) / 3, bw = bw,
                          from = lo, to = up, n = n)
      y <- 3 * d$y
    }
    data.frame(iter = g$iter[1], parameter = g$parameter[1],
               x = 10^d$x, density = y, ndensity = y / max(y))
  })
  do.call(rbind, out)
}

# true values as a data frame for geom_vline
abc_truth_to_df <- function(truth) {
  if (is.null(truth)) return(NULL)
  if (length(truth) != 5 || !is.numeric(truth)) {
    stop("truth should be 5 numbers: lambda0, lambda1, mu0, mu1, compl_rate")
  }
  if (is.null(names(truth))) names(truth) <- abc_param_names
  if (!setequal(names(truth), abc_param_names)) {
    # PBD::pbd_sim / pbd_sim_rcpp use different names (and mu1 means the
    # extinction rate of GOOD species there), so do not guess
    stop("truth should be named lambda0, lambda1, mu0, mu1, compl_rate ",
         "(the names used in the ABC output), or be unnamed in that order")
  }
  data.frame(parameter = factor(abc_param_labels[names(truth)],
                                levels = abc_param_labels),
             value = unname(truth))
}

# invisible points at the prior bounds, so each panel spans the prior
abc_prior_range_df <- function(lower, upper) {
  if (is.null(lower) || is.null(upper)) return(NULL)
  if (length(lower) != 5 || length(upper) != 5) {
    stop("lower and upper should each have 5 values (log10 prior bounds)")
  }
  data.frame(parameter = factor(rep(abc_param_labels, 2),
                                levels = abc_param_labels),
             value = 10^c(lower, upper))
}

#' Plot the posterior densities of the five PBD parameters
#'
#' Weighted density of each parameter in one generation of an ABC-SMC run,
#' on a log10 axis, with the true value as a dashed line if given.
#'
#' @param res output of \code{\link{perform_abc_par}}
#' @param truth optional true parameter values: 5 numbers, either named
#' lambda0, lambda1, mu0, mu1, compl_rate (the names in the ABC output) or
#' unnamed in that order. Note that \code{pbd_sim_rcpp} uses a different order
#' (b1, la1, b2, mu1, mu2), so convert first: lambda0 = b1, lambda1 = b2,
#' mu0 = mu1, mu1 = mu2, compl_rate = la1.
#' @param iteration generation to plot; defaults to the last one
#' @param lower,upper optional log10 prior bounds (as given to
#' perform_abc_par); if given, each panel spans the full prior range, which
#' shows how much the posterior has narrowed, and the densities are
#' restricted to the prior range (with reflection at the bounds)
#' @return a ggplot object
#' @examples
#' \dontrun{
#' plot_posterior(res, truth = c(lambda0 = 0.5, lambda1 = 0.5, mu0 = 0.1,
#'                               mu1 = 0.1, compl_rate = 0.3))
#' }
#' @export
plot_posterior <- function(res, truth = NULL, iteration = NULL,
                           lower = NULL, upper = NULL) {
  if (is.null(iteration)) iteration <- max(as.data.frame(res)$iter)
  long <- abc_to_long(res, iteration)
  truth_df <- abc_truth_to_df(truth)
  range_df <- abc_prior_range_df(lower, upper)

  dens <- abc_weighted_density(long, lower, upper)

  p <- ggplot2::ggplot(dens, ggplot2::aes(x = x)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = 0, ymax = density),
                         fill = "#5DCAA5", alpha = 0.6) +
    ggplot2::geom_line(ggplot2::aes(y = density), colour = "#0F6E56") +
    ggplot2::facet_wrap(~parameter, scales = "free", ncol = 3) +
    ggplot2::scale_x_log10() +
    ggplot2::labs(x = "parameter value (log scale)",
                  y = "posterior density (per log10 unit)",
                  title = paste("ABC posterior, generation", iteration)) +
    ggplot2::theme_bw()
  if (!is.null(range_df)) {
    p <- p + ggplot2::geom_blank(data = range_df,
                                 ggplot2::aes(x = value), inherit.aes = FALSE)
  }
  if (!is.null(truth_df)) {
    p <- p + ggplot2::geom_vline(data = truth_df,
                                 ggplot2::aes(xintercept = value),
                                 colour = "#D85A30", linetype = "dashed",
                                 linewidth = 0.8)
  }
  p
}

#' Plot how the parameter densities change over the ABC generations
#'
#' Stacked (ridgeline) weighted densities of each parameter, one ridge per
#' generation, from the first generation at the top to the last at the
#' bottom, on a log10 axis, with the true value as a dashed line if given.
#' Each ridge is scaled to the same maximum height.
#'
#' @inheritParams plot_posterior
#' @param iterations generations to plot; defaults to all
#' @param palette viridis colour map for the generations (light = early, dark
#' = late): "mako" (default), "rocket", "viridis", "magma", "plasma",
#' "inferno", "cividis" or "turbo"
#' @return a ggplot object
#' @examples
#' \dontrun{
#' plot_progression(res, truth = c(lambda0 = 0.5, lambda1 = 0.5, mu0 = 0.1,
#'                                 mu1 = 0.1, compl_rate = 0.3))
#' }
#' @export
plot_progression <- function(res, truth = NULL, iterations = NULL,
                             lower = NULL, upper = NULL, palette = "mako") {
  long <- abc_to_long(res, iterations)
  truth_df <- abc_truth_to_df(truth)
  range_df <- abc_prior_range_df(lower, upper)

  # ggridges' own density stat ignores weights, so the weighted densities
  # are computed here and drawn as ridgelines (each scaled to maximum 1)
  dens <- abc_weighted_density(long, lower, upper)

  # first generation at the top, last at the bottom
  dens$generation <- factor(dens$iter, levels = rev(sort(unique(dens$iter))))

  p <- ggplot2::ggplot(dens,
                       ggplot2::aes(x = x, y = generation, height = ndensity,
                                    fill = iter, group = iter)) +
    ggridges::geom_ridgeline(scale = 0.95, alpha = 0.8, colour = "grey30",
                             linewidth = 0.3) +
    ggplot2::facet_wrap(~parameter, scales = "free_x", nrow = 1) +
    ggplot2::scale_x_log10() +
    # light for early generations, dark for late ones; the extremes of the
    # map are skipped so no ridge is near-white or near-black
    ggplot2::scale_fill_viridis_c(option = palette, begin = 0.2, end = 0.85,
                                  direction = -1, guide = "none") +
    ggplot2::labs(x = "parameter value (log scale)", y = "generation",
                  title = "ABC posterior over the generations") +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.spacing.x = ggplot2::unit(1.2, "lines"))
  if (!is.null(range_df)) {
    p <- p + ggplot2::geom_blank(data = range_df, ggplot2::aes(x = value),
                                 inherit.aes = FALSE)
  }
  if (!is.null(truth_df)) {
    p <- p + ggplot2::geom_vline(data = truth_df,
                                 ggplot2::aes(xintercept = value),
                                 colour = "#D85A30", linetype = "dashed",
                                 linewidth = 0.8)
  }
  p
}
