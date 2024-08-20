set.seed(42)
#set.seed(5)
#set.seed(3)
test_tree <- ape::rbdtree(birth = 1, death = 0.0, Tmax = 5)

g_tree <- treestats::gamma_statistic(test_tree)
col_tree <- treestats::colless(test_tree)
num_lin <- treestats::number_of_lineages(test_tree)

cat(g_tree, col_tree, num_lin, "\n")

res <- pbdABC::test_abc_par(num_particles = 1000,
                            num_iterations = 15,
                            crown_age = 5,
                            min_lin = num_lin*0.5,
                            max_lin = num_lin*2,
                            lambdas = c(1, 10, 10000, 10000, 0.01),
                            obs_gamma = g_tree,
                            obs_colless = col_tree,
                            obs_num_lin = num_lin,
                            num_threads = 1)
