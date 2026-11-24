#pragma once

#include "prior.h"
#include "rnd_thijs.h"
#include "sim_pbd.h"
#include "ltable.h"

#include "L2newick.h"

struct particle {
  param_set params_;
  double gamma;
  double colless;
  int num_lin;

  bool success = false;

  double weight = 1.0;
  double prob_prior;

  std::string newick_tree;

  particle(prior prior_dist, rnd_t& rndgen) {
    params_ = prior_dist.gen_prior(rndgen);
   // params_.back() = 1.0 / params_.back();
    success = false;
  }

  particle() {
    params_ = {0.5, 0.1, 0.5, 0.1, 1e6};
    success = false;
  }

  particle(const particle& other) = default;
  particle& operator=(const particle& other) = default;

  // Perturb all parameters jointly with independent normal steps on the log
  // scale, so the transition density is the product computed in prob_perturb.
  void perturb(rnd_t& rndgen) {
    for (size_t i = 0; i < params_.size(); ++i) {
      params_[i] = rndgen.perturb_particle_val(params_[i], i);
    }
    return;
  }

  double prob_perturb(const particle& other,
                      const rnd_t& rndgen) {

    // Density of moving from 'other' to this particle in log-parameter space:
    // product of independent normal densities, accumulated as a log sum.
    const double log_2pi = std::log(2.0 * 3.141592653589793238);
    double log_k = 0.0;
    for (size_t i = 0; i < other.params_.size(); ++i) {
      double sigma = rndgen.kernel_sigmas[i];
      double d = std::log(params_[i]) - std::log(other.params_[i]);
      log_k += -0.5 * (d * d) / (sigma * sigma) - std::log(sigma) - 0.5 * log_2pi;
    }
    return std::exp(log_k);
  }

  bool pass_prior(const prior& prior_dist) {
      auto param_check = params_;
     // param_check.back() = 1.0 / param_check.back();
      prob_prior = prior_dist.dens_prior(param_check);
      if (prob_prior > 0.0) return true;

      return false;
  }

  void update_weight(const std::vector<particle>& other,
                     const rnd_t& rndgen,
                     const prior& prior_dist) {
    double sum_perturb = 0.0;
    for (const auto& i : other) {
      double prob = prob_perturb(i, rndgen);
      if (std::isnan(prob)) prob = 0.0; // this can happen sometimes for impossible values
      if (std::isinf(prob)) prob = 0.0;
      sum_perturb += prob * i.weight;
    }

    auto prior_check = params_;
   // prior_check.back() = 1.0 / prior_check.back();
    prob_prior = prior_dist.dens_prior(prior_check);

    // The kernel is a density over log(theta), so the prior must be too:
    // f_log(log theta) = f(theta) * prod(theta).
    double log_jacobian = 0.0;
    for (const auto& v : params_) log_jacobian += std::log(v);

    auto new_weight = prob_prior * std::exp(log_jacobian) / sum_perturb;
    weight = new_weight;
    if (std::isnan(new_weight)) weight = 0.0;
    if (std::isinf(new_weight)) weight = 1e10;
    if (sum_perturb == 0.0) weight = 0.0;
    if (prob_prior == 0.0) weight = 0.0;
  }

  void sim(double crown_age,
           int min_lin,
           int max_lin) {

    sim_pbd sim(params_, crown_age, max_lin * 10);
    sim.run();

    colless = 1e6; // default bad values
    gamma = 1e6;
    success = false;
    num_lin = 1e6;

    if (sim.status != "success") return;

    // the number of species in the reconstructed tree (= number of tips, one
    // sampled lineage per extant species), which is what an observed tree
    // gives; the number of good lineages is not the same thing. After
    // select_random, the sampled tips are exactly the rows still alive.
    int num_tips = 0;
    for (const auto& row : sim.L) {
      if (row[species_property::death_time] < 0) num_tips++;
    }
    num_lin = num_tips;

    if (num_lin    >= min_lin   &&
        num_lin    <= max_lin   &&
        num_lin    >= 2) {
        auto ltable_ = drop_extinct(sim.L);
        bool crowns_alive = ltable_[0][species_property::death_time] < 0 &&
                            ltable_[1][species_property::death_time] < 0;

        double resulting_crown_age = ltable_[0][species_property::birth_time];

        if (crowns_alive == true &&
            resulting_crown_age == crown_age)  {
          std::vector<double> brts = brts_from_ltable(ltable_);
          gamma = calc_gamma(brts);

          colless_stat_ltable s(ltable_);
          colless = static_cast<double>(s.colless());
          success = true;

          newick_tree = ltable_to_newick(ltable_,
                                         true);
        } else {
          num_lin = 1e6; // to trigger bad fit
        }
    } else {
      num_lin = 1e6; // to trigger bad fit
    }
    return;
  }
};

// Beaumont et al. (2009): kernel variance = 2 * weighted variance of the
// current population (here on the log scale). Returns standard deviations.
inline std::array<double, 5> compute_kernel_sigmas(const std::vector<particle>& sample) {
  std::array<double, 5> sigmas = {1.0, 1.0, 1.0, 1.0, 1.0};
  if (sample.empty()) return sigmas;

  double sum_w = 0.0;
  for (const auto& p : sample) sum_w += p.weight;
  const bool use_weights = sum_w > 0.0 && std::isfinite(sum_w);

  for (size_t j = 0; j < 5; ++j) {
    double mean = 0.0;
    for (const auto& p : sample) {
      double w = use_weights ? p.weight / sum_w : 1.0 / sample.size();
      mean += w * std::log(p.params_[j]);
    }
    double var = 0.0;
    for (const auto& p : sample) {
      double w = use_weights ? p.weight / sum_w : 1.0 / sample.size();
      double d = std::log(p.params_[j]) - mean;
      var += w * d * d;
    }
    double s = std::sqrt(2.0 * var);
    if (!std::isfinite(s) || s < 1e-6) s = 1e-6;  // avoid a degenerate kernel
    sigmas[j] = s;
  }
  return sigmas;
}
