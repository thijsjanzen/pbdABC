#pragma once

#include <vector>
#include <array>
#include <cmath>

#include "prior.h"
#include "rnd_thijs.h"
#include "sim_pbd.h"
#include "ltable.h"
#include "gamma.h"
#include "colless.h"
#include "particle.h"

struct analysis {
  std::vector<particle> current_sample;

  std::vector<particle> new_sample;

  rnd_t rndgen_;

  const double ref_gamma;
  const double ref_colless;
  const double ref_num_lin;
  const double crown_age;
  const double min_lin;
  const double max_lin;
  const int num_particles;

  prior prior_dist;

  std::vector<double> threshold;

  analysis(int n,
           int num_iterations,
           double ca,
           double minimum_lineages,
           double maximum_lineages,
           param_set lower,
           param_set upper,
           param_set means,
           double obs_gamma,
           double obs_colless,
           double obs_num_lin) :
    ref_gamma(obs_gamma),
    ref_colless(obs_colless),
    ref_num_lin(obs_num_lin),
    crown_age(ca),
    min_lin(minimum_lineages),
    max_lin(maximum_lineages),
    num_particles(n) {
    rndgen_ = rnd_t();
    for (size_t i = 0; i < num_iterations; ++i) {
      threshold.push_back(10 * std::exp(-0.5 * (i - 1)));
    }

    if (means.empty()) {
      prior_dist = prior(lower, upper);
    } else if (means[0] < 0) {
      prior_dist = prior(lower, upper);
    } else {
      prior_dist = prior(means);
    }
  }

  void iterate_first() {

    std::cerr << "0--------25--------50--------75--------100\n";
    std::cerr << "*";
    int updateFreq = num_particles / 20;
    if(updateFreq < 1) updateFreq = 1;

    int prev_print = 0;
    while(current_sample.size() < num_particles) {
      auto new_particle = particle(prior_dist, rndgen_);

      new_particle.sim(crown_age, min_lin, max_lin);

      if (new_particle.num_lin >= min_lin &&
          new_particle.num_lin <= max_lin) {
        new_particle.success = true;
        current_sample.push_back(new_particle);
      }

      if (current_sample.size() % updateFreq == 0 &&
          current_sample.size() != prev_print) {
        std::cerr << "**";
        prev_print = current_sample.size();
      }
    }

    std::cerr << "\n";
  }

  void iterate(int iteration) {
    new_sample.clear();
    std::vector<double> weights (current_sample.size());
    for (size_t i = 0; i < current_sample.size(); ++i) {
      weights[i] = current_sample[i].weight;
    }

    std::discrete_distribution pick_particle(weights.begin(), weights.end());

    std::cerr << "0--------25--------50--------75--------100\n";
    std::cerr << "*";
    int updateFreq = num_particles / 20;
    if(updateFreq < 1) updateFreq = 1;
    size_t prev_print = 0;

    while(new_sample.size() < num_particles) {
        auto new_particle = current_sample[pick_particle(rndgen_.rndgen_)];
        new_particle.perturb(rndgen_);

        if (new_particle.pass_prior(prior_dist)) {
          new_particle.sim(crown_age, min_lin, max_lin);
          if (new_particle.success == true) {

            double dist = calc_dist(new_particle);
            if (dist < threshold[iteration]) {
              new_particle.update_weight(current_sample, rndgen_, prior_dist);
              new_sample.push_back(new_particle);
            }
          }
        }

        if (new_sample.size() % updateFreq == 0 &&
            new_sample.size() != prev_print) {
          std::cerr << "**";
          prev_print = new_sample.size();
        }
    }

    current_sample = new_sample;
    std::cerr << "\n";
  }


  double calc_dist(const particle& p) {
    double d1 = (p.gamma - ref_gamma);
    double d2 = (p.colless - ref_colless);
    double d3 = (p.num_lin- ref_num_lin);

    std::array<double, 3> diff;
    diff[0] = (d1 * d1) / ref_gamma;
    diff[1] = (d2 * d2) / ref_colless;
    diff[2] = (d3 * d3) * 1.0 / ref_num_lin;

    return *std::max_element(diff.begin(), diff.end());
  }
};
