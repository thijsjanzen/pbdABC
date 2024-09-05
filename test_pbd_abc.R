lambda0 <- 0.2
compl_rate <- 0.001
lambda1 <- 0.2
mu0 <- 0
mu1 <- 0

crown_age <- 48

#found <- c()
#for (i in 1:100) {
#  test_tree <- pbdABC::pbd_sim_rcpp(pars = c(lambda0, compl_rate, lambda1, mu0, mu1),
#                                  age = crown_age)
#  found[i] <- treestats::number_of_lineages(test_tree)
#}
#median(found)
#hist(found)
#mean(found)

test_tree <- pbdABC::pbd_sim_rcpp(pars = c(lambda0, compl_rate, lambda1, mu0, mu1),
                                  age = crown_age)
treestats::number_of_lineages(test_tree)

g_tree <- treestats::gamma_statistic(test_tree)
col_tree <- treestats::colless(test_tree)
num_lin <- treestats::number_of_lineages(test_tree)

cat(g_tree, col_tree, num_lin, "\n")

estim <- PBD::pbd_ML(treestats::branching_times(test_tree),
                     initparsopt = c(lambda0, compl_rate),
                     idparsopt = c(1, 3),
                     idparsfix = c(2, 4),
                     parsfix = c(0, 0),
                     exteq = 0,
                     cond = 1,
                     btorph = 1,
                     soc = 2)

res <- pbdABC::perform_abc_par(num_particles = 1000,
                               num_iterations = 15,
                               crown_age = crown_age,
                               min_lin = num_lin * 0.5,
                               max_lin = num_lin * 2,
                               lower = c(-2, -2, -6, -6, -4),
                               upper = c(0, 0, -5, -5, 1),
                               #means = c(1 / lambda0, 1/lambda1, 1000, 1000, 1 / compl_rate),
                               obs_gamma = g_tree,
                               obs_colless = col_tree,
                               obs_num_lin = num_lin,
                               num_threads = 10,
                               sim_number = 1,
                               bd_lambda = estim$b,
                               bd_mu = 0.01)










if (1 == 2) {


res <- pbdABC::perform_abc_par(num_particles = 100,
                               num_iterations = 15,
                               crown_age = crown_age,
                               min_lin = num_lin * 0.5,
                               max_lin = num_lin * 2,
                               means = c(1 / lambda0, 1/lambda1, 1000, 1000, 1 / compl_rate),
                               obs_gamma = g_tree,
                               obs_colless = col_tree,
                               obs_num_lin = num_lin,
                               num_threads = 10,
                               sim_number = 1,
                               bd_lambda = estim$b,
                               bd_mu = 0.01)
}




