#pragma once

#include "prior.h"
#include "rnd_thijs.h"
#include "sim_pbd.h"
#include "ltable.h"

struct particle {
  param_set params_;
  double gamma;
  double colless;
  int num_lin;

  bool success = false;

  double weight = 1.0;
  ltable ltable_;
  double prob_prior;

  particle(prior prior_dist, rnd_t& rndgen) {
    params_ = prior_dist.gen_prior(rndgen);
    success = false;
  }

  particle() {
    params_ = {0.5, 0.1, 0.5, 0.1, 1e6};
    success = false;
  }

  particle(const particle& other) = default;
  particle& operator=(const particle& other) = default;

  void perturb(rnd_t& rndgen) {
    size_t index = rndgen.random_number(params_.size());
    params_[index] = rndgen.perturb_particle_val(params_[index], index);
    return;
  }

  double prob_perturb(const particle& other,
                      const rnd_t& rndgen) {

    // alternative calculation
    double alt_prob = 0.0;
    for (size_t i = 0; i < other.params_.size(); ++i) {
      double sigma = rndgen.kernel_sigmas[i];
      double factor = 1.0 / (-2 * sigma * sigma);
      double prefactor = 1.0 / std::sqrt(2 * 3.141592653589793238 * sigma * sigma);

      double d = std::log(params_[i]) - std::log(other.params_[i]);
      double p = (d * d) * factor;
      alt_prob += std::exp(p) * prefactor;
    }
    return alt_prob;
  }

  bool pass_prior(const prior& prior_dist) {
      prob_prior = prior_dist.dens_prior(params_);
      if (prob_prior > 0.0) return true;

      return false;
  }

  void update_weight(const std::vector<particle>& other,
                     const rnd_t& rndgen,
                     const prior& prior_dist) {
    double sum_perturb = 0.0;
    for (const auto& i : other) {
      double prob = prob_perturb(i, rndgen);
      sum_perturb += prob * i.weight;
    }

    prob_prior = prior_dist.dens_prior(params_);
    auto new_weight = prob_prior / sum_perturb;
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
    num_lin = sim.num_good_species;

    colless = 1e6; // default bad values
    gamma = 1e6;
    success = false;

    if (sim.status == "success" &&
        num_lin    >= min_lin   &&
        num_lin    <= max_lin) {

        ltable_ = drop_extinct(sim.L);
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
        } else {
          num_lin = 1e6; // to trigger bad fit
        }
    } else {
      num_lin = 1e6; // to trigger bad fit
    }
    return;
  }
};
