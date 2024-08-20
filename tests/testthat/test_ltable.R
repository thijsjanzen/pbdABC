context("ltable")

test_that("usage", {

  test_result <- function(focal_tree) {
    tree1 <- treestats::l_to_phylo(focal_tree$ltable, drop_extinct = TRUE)
    tree2 <- treestats::l_to_phylo(focal_tree$ltable_dropped, drop_extinct = TRUE)
    s1 <- treestats::calc_all_stats(tree1)
    s2 <- treestats::calc_all_stats(tree2)
    testthat::expect_equal(s1, s2)

    g1 <- treestats::gamma_statistic(tree1)
    g2 <- treestats::gamma_statistic(tree2)
    g3 <- treestats::gamma_statistic(focal_tree$ltable_dropped)

    testthat::expect_equal(g1, g2)
    testthat::expect_equal(g1, g3)
    testthat::expect_equal(g1, focal_tree$gamma)

    c1 <- treestats::colless(tree1)
    c2 <- treestats::colless(tree2)
    c3 <- treestats::colless(focal_tree$ltable_dropped)

    testthat::expect_equal(c1, c2)
    testthat::expect_equal(c1, c3)
    testthat::expect_equal(c1, focal_tree$colless)
  }

  for (r in 1:10) {
    local_tree <- pbdABC::test_simulations(birth = 1, death = 0.0, crown_age = 2)
    test_result(local_tree)

    local_tree <- pbdABC::test_simulations(birth = 1, death = 0.1, crown_age = 5)
    test_result(local_tree)

    local_tree <- pbdABC::test_simulations(birth = 1, death = 0.4, crown_age = 5)
    test_result(local_tree)
  }
})
