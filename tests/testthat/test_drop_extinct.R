context("drop_extinct")

test_that("dropped L table matches treestats on protracted trees", {
  # the existing ltable tests use (near-)instant speciation completion; here
  # completion is slow, so many extant lineages are not the sampled
  # representative of their species and get an artificial death, including,
  # regularly, the crown lineages themselves
  set.seed(42)
  pars <- list(c(0.6, 0.1, 0.6, 0.1, 0.5),   # la0, mu0, la1, mu1, compl
               c(1.0, 0.8, 0.5, 0.4, 1.0),
               c(0.5, 0.4, 2.0, 0.3, 0.1))
  num_checked <- 0
  for (p in pars) {
    for (r in 1:30) {
      x <- pbdABC:::sim_pbd_cpp(p[1], p[2], p[3], p[4], p[5], 5, 2000, 1)
      if (x$status != "success") {
        testthat::expect_equal(nrow(x$ltable_dropped), 0)
        next
      }
      if (sum(x$ltable[, 4] == -1) < 3) next
      ref <- treestats::l_to_phylo(x$ltable, drop_extinct = TRUE)
      got <- treestats::l_to_phylo(x$ltable_dropped, drop_extinct = TRUE)
      testthat::expect_true(ape::all.equal.phylo(ref, got,
                                                 use.tip.label = FALSE))
      testthat::expect_equal(treestats::gamma_statistic(x$ltable_dropped),
                             treestats::gamma_statistic(ref))
      testthat::expect_equal(treestats::colless(x$ltable_dropped),
                             treestats::colless(ref))
      num_checked <- num_checked + 1
    }
  }
  testthat::expect_gt(num_checked, 20)
})
