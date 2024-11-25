#include <Rcpp.h>
#include <fstream>
#include "abc.h"
#include "util.h"

#include "abc_par.h"

size_t get_rcpp_num_threads_abc() {
  auto* nt_env = std::getenv("RCPP_PARALLEL_NUM_THREADS");
  return (nullptr == nt_env)
    ? tbb::task_arena::automatic  // -1
  : static_cast<size_t>(std::atoi(nt_env));
}


//' function to do abc using rcpp
//' @param num_particles number of particles
//' @param num_iterations number of iterations
//' @param crown_age crown age
//' @param min_lin minimum number of lineages from the prior
//' @param max_lin maximum number oflineages from the prior
//' @param lower minimum values of log-uniform prior (e.g. -3 = 10^-3), ordering:
//' l0, l1, mu0, mu1, compl_rate
//' @param use_inv_prior usage of an inverse prior on completion rate?
//' @param upper upper values of log-uniform prior (e.g. 3 = 10^3)//' @param obs_gamma observed gamma value to fit on
//' @param obs_colless observed colless value to fit on
//' @param obs_num_lin observed number of lineages to fit on
//' @param sim_number internal value for Thijs
//' @param bd_lambda bd lambda
//' @param bd_mu bd mu
//' @export
//' @description
//' Fit to the data is assessed by the sum of differences for
//' [gamma, colless, num_lineages], rescaled by the observed value, e.g.
//' ((O-E)^2)/E, where O is the value of the proposed simulation and E is the
//' value of the empirical data (e.g. obs_gamma, obs_colless or obs_num_lineages).
//' The acceptance threshold diminishes exponentially.
// [[Rcpp::export]]
Rcpp::NumericMatrix perform_abc_rcpp_par(int num_particles,
                                         int num_iterations,
                                         double crown_age,
                                         double min_lin,
                                         double max_lin,
                                         std::vector<double> lower,
                                         std::vector<double> upper,
                                         std::vector<double> means,
                                         bool use_inv_prior,
                                         double obs_gamma,
                                         double obs_colless,
                                         double obs_num_lin,
                                         int sim_number,
                                         double bd_lambda,
                                         double bd_mu,
                                         double limiting_accept_rate) {

   auto num_threads = get_rcpp_num_threads_abc();
   auto global_control = tbb::global_control(tbb::global_control::max_allowed_parallelism, num_threads);

   std::array<double, 5> low; std::copy(lower.begin(), lower.end(), low.begin());
   std::array<double, 5> up;  std::copy(upper.begin(), upper.end(), up.begin());
   std::array<double, 5> m;   std::copy(means.begin(), means.end(), m.begin());

   analysis_par focal_analysis(num_particles,
                               num_iterations,
                               crown_age,
                               min_lin,
                               max_lin,
                               low,
                               up,
                               m,
                               use_inv_prior,
                               obs_gamma,
                               obs_colless,
                               obs_num_lin,
                               limiting_accept_rate);

   // check if is first!

   focal_analysis.find_output(sim_number);
   std::vector< std::array<double, 10>> res;

   if (focal_analysis.blank_sheet) {
     focal_analysis.iterate_first();
     update_output(res, focal_analysis.current_sample, 0);
     focal_analysis.write_to_file(sim_number, 0,
                                  obs_gamma, obs_colless, obs_num_lin,
                                  bd_lambda, bd_mu);
     focal_analysis.write_trees(sim_number, 0);
     focal_analysis.update_kernel(0);
   }

   size_t i = 1;
   if (!focal_analysis.blank_sheet) i = focal_analysis.starting_iteration;

   for (; i < num_iterations; ++i) {
     focal_analysis.iterate(i);

     update_output(res, focal_analysis.current_sample, i);

     focal_analysis.write_to_file(sim_number, i,
                                  obs_gamma, obs_colless, obs_num_lin,
                                  bd_lambda, bd_mu);
     focal_analysis.write_trees(sim_number, i);
     focal_analysis.update_kernel(i);

     if (focal_analysis.accept_rate < limiting_accept_rate) break;
   }

   Rcpp::NumericMatrix out;
   particle_to_numericmatrix(res, out);

   return out;
 }

//' function to do abc using rcpp
//' @param num_particles number of particles
//' @param num_iterations number of iterations
//' @param crown_age crown age
//' @param min_lin minimum number of lineages from the prior
//' @param max_lin maximum number oflineages from the prior
//' @param lower minimum values of log-uniform prior (e.g. -3 = 10^-3), ordering:
//' l0, l1, mu0, mu1, compl_rate
//' @param upper upper values of log-uniform prior (e.g. 3 = 10^3)//' @param obs_gamma observed gamma value to fit on
//' @param use_inv_prior use inverse prior
//' @param obs_gamma observed gamma value to fit on
//' @param obs_colless observed colless value to fit on
//' @param obs_num_lin observed number of lineages to fit on
//' @description
//' Fit to the data is assessed by the sum of differences for
//' [gamma, colless, num_lineages], rescaled by the observed value, e.g.
//' ((O-E)^2)/E, where O is the value of the proposed simulation and E is the
//' value of the empirical data (e.g. obs_gamma, obs_colless or obs_num_lineages).
//' The acceptance threshold diminishes exponentially.
// [[Rcpp::export]]
double test_abc_rcpp_par(int num_particles,
                                      int num_iterations,
                                      double crown_age,
                                      double min_lin,
                                      double max_lin,
                                      std::vector<double> lower,
                                      std::vector<double> upper,
                                      std::vector<double> means,
                                      bool use_inv_prior,
                                      double obs_gamma,
                                      double obs_colless,
                                      double obs_num_lin) {

   auto clock_start = std::chrono::system_clock::now();


   auto num_threads = get_rcpp_num_threads_abc();
   auto global_control = tbb::global_control(tbb::global_control::max_allowed_parallelism, num_threads);

   std::array<double, 5> low; std::copy(lower.begin(), lower.end(), low.begin());
   std::array<double, 5> up;  std::copy(upper.begin(), upper.end(), up.begin());
   std::array<double, 5> m;   std::copy(means.begin(), means.end(), m.begin());

   analysis_par focal_analysis(num_particles,
                               num_iterations,
                               crown_age,
                               min_lin,
                               max_lin,
                               low,
                               up,
                               m,
                               use_inv_prior,
                               obs_gamma,
                               obs_colless,
                               obs_num_lin,
                               1e-6);

   std::vector< std::array<double, 10>> res;

   focal_analysis.iterate_first();
   update_output(res, focal_analysis.current_sample, 0);
   focal_analysis.update_kernel(0);

   focal_analysis.iterate(1);
   update_output(res, focal_analysis.current_sample, 1);
   focal_analysis.update_kernel(1);



   Rcpp::NumericMatrix out;
   particle_to_numericmatrix(res, out);

   auto clock_now = std::chrono::system_clock::now();
   std::chrono::duration<double> elapsed_seconds = clock_now - clock_start;
   std::cout << "this took: " << elapsed_seconds.count() << "seconds\n";


   return elapsed_seconds.count();
}




//' function to do abc using rcpp
//' @param num_particles number of particles
//' @param num_iterations number of iterations
//' @param crown_age crown age
//' @param min_lin minimum number of lineages from the prior
//' @param max_lin maximum number oflineages from the prior
//' @param lower minimum values of log-uniform prior (e.g. -3 = 10^-3), ordering:
//' l0, l1, mu0, mu1, compl_rate
//' @param upper upper values of log-uniform prior (e.g. 3 = 10^3)
//' @param obs_gamma observed gamma value to fit on
//' @param obs_gamma observed gamma value to fit on
//' @param obs_colless observed colless value to fit on
//' @param obs_num_lin observed number of lineages to fit on
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix perform_abc_rcpp(int num_particles,
                                      int num_iterations,
                                      double crown_age,
                                      double min_lin,
                                      double max_lin,
                                      std::vector<double> lower,
                                      std::vector<double> upper,
                                      std::vector<double> means,
                                      double obs_gamma,
                                      double obs_colless,
                                      double obs_num_lin) {

   std::array<double, 5> low; std::copy(lower.begin(), lower.end(), low.begin());
   std::array<double, 5> up;  std::copy(upper.begin(), upper.end(), up.begin());
   std::array<double, 5> m;   std::copy(means.begin(), means.end(), m.begin());

   analysis focal_analysis(num_particles,
                           num_iterations,
                           crown_age,
                           min_lin,
                           max_lin,
                           low,
                           up,
                           m,
                           obs_gamma,
                           obs_colless,
                           obs_num_lin);

   std::vector< std::array<double, 10>> res;

   focal_analysis.iterate_first();
   update_output(res, focal_analysis.current_sample, 0);

   for (size_t i = 1; i < num_iterations; ++i) {
     focal_analysis.iterate(i);
     update_output(res, focal_analysis.current_sample, i);
   }

   Rcpp::NumericMatrix out;
   particle_to_numericmatrix(res, out);

   return out;
 }



//' function to test prior distribution generation
//' @param birth birth
//' @param death death
//' @param crown_age crown age
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix test_prior(int num_samples,
                               const std::vector<double>& means,
                               bool use_inv) {

   std::vector< std::array< double, 2> > generated(num_samples);

   rnd_t rndgen;
   std::array<double, 5> m;   std::copy(means.begin(), means.end(), m.begin());

   prior prior_dist(m, use_inv);

   for (int i = 0; i < num_samples;  ++i) {
      auto val = prior_dist.gen_prior(rndgen);
      generated[i][0] = val.back();

      if (use_inv) {
         generated[i][1] = prior_dist.dens_inv_only_compl_exp(val);
      } else {
         generated[i][1] = prior_dist.dens_only_compl_exp(val);
      }
   }

   Rcpp::NumericMatrix out(num_samples, 2);

   for (int i = 0; i < num_samples; ++i) {
      for (int j = 0; j < 2; ++j) {
         out(i, j) = generated[i][j];
      }
   }

   return out;
}




//' function to test drop extinct and calculate stats
//' @param birth birth
//' @param death death
//' @param crown_age crown age
//' @export
// [[Rcpp::export]]
Rcpp::List test_simulations(double birth,
                             double death,
                             double crown_age) {

   sim_pbd sim({birth, death, 0.0, 0.0, 1e6}, // l0, mu0, l1, mu1, tau
               crown_age, 300 * 10);

   bool keep_running = true;
   int num_tries = 0;
   while(keep_running) {

     sim.run();

     int num_lin = 0;
     for (const auto& i : sim.L) {
       if (i[species_property::death_time] < 0) num_lin++;
     }
     bool crowns_alive = sim.L[0][species_property::death_time] < 0 &&
       sim.L[1][species_property::death_time] < 0;

     if (crowns_alive && num_lin > 10) {
       keep_running = false;
     }
     num_tries++;
     if (num_tries > 10000) {
       Rcpp::stop("Can't complete simulation");
     }
   }

   auto ltable_ = drop_extinct(sim.L);

   std::vector<double> brts = brts_from_ltable(ltable_);

   auto gamma = calc_gamma(brts);

   colless_stat_ltable s(ltable_);
   auto colless = static_cast<double>(s.colless());

   Rcpp::NumericMatrix raw_ltable;
   vector_to_numericmatrix(sim.L, raw_ltable);

   Rcpp::NumericMatrix dropped_ltable;
   vector_to_numericmatrix(ltable_, dropped_ltable);

   return Rcpp::List::create(Rcpp::Named("ltable") = raw_ltable,
                             Rcpp::Named("ltable_dropped") = dropped_ltable,
                             Rcpp::Named("gamma") = gamma,
                             Rcpp::Named("colless") = colless);
 }
