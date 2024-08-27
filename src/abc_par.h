#pragma once

#include <vector>
#include <array>
#include "colless.h"
#include "gamma.h"

#include "rnd_thijs.h"
#include "sim_pbd.h"
#include "ltable.h"

#include <mutex>

#include <RcppParallel.h>

struct particle_par {
  std::array<double, 5> params_;
  double gamma;
  double colless;
  int num_lin;

  bool success = false;

  double weight = 1.0;
  ltable ltable_;

  particle_par(rnd_t& rndgen) {
    params_ = rndgen.draw_from_prior();
    success = false;
  }

  particle_par() {
    params_ = {0.5, 0.1, 0.5, 0.1, 1e6};
    success = false;
  }

  particle_par(const particle_par& other) = default;
  particle_par& operator=(const particle_par& other) = default;

  void perturb(rnd_t& rndgen) {
    size_t index = rndgen.random_number(params_.size());
    params_[index] = rndgen.perturb_particle_val(params_[index], index);
    return;
  }

  double prob_perturb(const particle_par& other,
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

  void update_weight(const std::vector<particle_par>& other,
                     const rnd_t& rndgen) {
    double sum_perturb = 0.0;
    for (const auto& i : other) {
      double prob = prob_perturb(i, rndgen);
      sum_perturb += prob * i.weight;
    }

    double prob_prior = rndgen.dens_prior(params_);
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
    num_lin = 0;
    for (const auto& i : sim.L) {
      if (i[3] < 0) num_lin++;
    }
    bool crowns_alive = sim.L[0][species_property::death_time] < 0 &&
      sim.L[1][species_property::death_time] < 0;

    if (sim.status == "success" &&
        num_lin >= min_lin &&
        num_lin <= max_lin &&
        crowns_alive) {

      ltable_ = drop_extinct(sim.L);

      std::vector<double> brts = brts_from_ltable(ltable_);
      gamma = calc_gamma(brts);

      colless_stat_ltable s(ltable_);
      colless = static_cast<double>(s.colless());
      success = true;
    } else {
      colless = 1e6;
      gamma = 1e6;
      num_lin = 1e6;
      success = false;
    }
    return;
  }
};

struct analysis_par {
  std::vector<particle_par> current_sample;

  std::vector<particle_par> new_sample;

  rnd_t rndgen_;

  const double ref_gamma;
  const double ref_colless;
  const double ref_num_lin;
  const double crown_age;
  const double min_lin;
  const double max_lin;
  const double limit_accept_rate;
  const int num_particles;
  double accept_rate;

  std::vector<double> threshold;
  std::array<double, 5> sigmas;

  analysis_par(int n,
           int num_iterations,
           double ca,
           double minimum_lineages,
           double maximum_lineages,
           std::vector<double> lower,
           std::vector<double> upper,
           double obs_gamma,
           double obs_colless,
           double obs_num_lin,
           double limit_rate) :
    ref_gamma(obs_gamma),
    ref_colless(obs_colless),
    ref_num_lin(obs_num_lin),
    crown_age(ca),
    min_lin(minimum_lineages),
    max_lin(maximum_lineages),
    num_particles(n),
    limit_accept_rate(limit_rate) {
    rndgen_ = rnd_t(lower, upper);
    for (size_t i = 0; i < num_iterations; ++i) {
      threshold.push_back(10 * std::exp(-0.5 * (i - 1)));
    }
  }

  std::vector<int> get_seeds(int num_seeds) {
    rnd_t rndgen;
    std::vector<int> seed_values;
    for (int i = 0; i < num_seeds; ++i) {
      seed_values.push_back(rndgen.random_number(INT_MAX)); // large value
    }
    return seed_values;
  }

  void iterate_first() {
    int num_remaining = num_particles - current_sample.size();
    accept_rate = 1.0;

    const int num_seeds = 20;

    Rcpp::Rcout << "0--------25--------50--------75--------100\n";
    Rcpp::Rcout << "*";
    int updateFreq = num_particles / 20;
    if(updateFreq < 1) updateFreq = 1;

    int prev_update = 0;

    while(current_sample.size() < num_particles) {
      for (size_t t = prev_update; t < current_sample.size(); ++t) {
        if (t % updateFreq == 0) {
          Rcpp::Rcout << "**";
        }
      }
      prev_update = current_sample.size();


      int loop_size = num_remaining * 1.0 / accept_rate;

      int seed_index = 0;
      std::mutex mutex;
      std::vector< int > seed_values = get_seeds(num_seeds);

      std::vector< particle_par > found_particles(loop_size);

      tbb::parallel_for(
        tbb::blocked_range<unsigned>(0, loop_size),
        [&](const tbb::blocked_range<unsigned>& r) {

          rnd_t rndgen2(seed_values[seed_index],
                                   rndgen_);
          {
            std::lock_guard<std::mutex> m(mutex);
            seed_index++;
            if (seed_index >= num_seeds) { // just in case.
              for (int i = 0; i < num_seeds; ++i) {
                seed_values[i] = rndgen2.random_number(INT_MAX);
              }
              seed_index = 0;
            }
          }

          for (unsigned i = r.begin(); i < r.end(); ++i) {

            auto new_particle = particle_par(rndgen2);

            new_particle.sim(crown_age, min_lin, max_lin);

            if (new_particle.num_lin >= min_lin &&
                new_particle.num_lin <= max_lin) {
              new_particle.success = true;
              found_particles[i] = new_particle;
            }
          }
        });

      int num_accepted = 0;
      for (const auto& i : found_particles) {
        if (i.success && current_sample.size() < num_particles) {
          current_sample.push_back(i);
          num_accepted++;
        }
      }
      accept_rate = 1.0 * (1 + num_accepted) / found_particles.size();
  //    Rcpp::Rcout << "current_accept_rate: " << accept_rate << "\n";
    }

    Rcpp::Rcout << "\n";
  }

  void iterate(int iteration) {
    new_sample.clear();
    std::vector<double> weights (current_sample.size());
    for (size_t i = 0; i < current_sample.size(); ++i) {
      weights[i] = current_sample[i].weight;
    }

    int num_remaining = num_particles;
    accept_rate = 1.0;
    const int num_seeds = 20;
    Rcpp::Rcout << "iteration: " << iteration << "\n";
    Rcpp::Rcout << "0--------25--------50--------75--------100\n";
    Rcpp::Rcout << "*";
    int updateFreq = num_particles / 20;
    if(updateFreq < 1) updateFreq = 1;

    int prev_update = 0;
    int num_tried = 0;
    while(new_sample.size() < num_particles) {
      for (size_t t = prev_update; t < new_sample.size(); ++t) {
        if (t % updateFreq == 0) {
          Rcpp::Rcout << "**";
        }
      }
      prev_update = new_sample.size();

      num_remaining = num_particles - new_sample.size();

      int loop_size = num_remaining * 1.0 / accept_rate;
      if (loop_size > 100000) loop_size = 100000;

      std::cerr << "\nloop_size: " << loop_size << "\n";


      int seed_index = 0;
      std::mutex mutex;
      std::vector< int > seed_values = get_seeds(num_seeds);

      std::vector< particle_par > found_particles(loop_size);

      tbb::parallel_for(
        tbb::blocked_range<unsigned>(0, loop_size),
        [&](const tbb::blocked_range<unsigned>& r) {

          rnd_t rndgen2(seed_values[seed_index],
                                   rndgen_);
          {
            std::lock_guard<std::mutex> m(mutex);
            seed_index++;
            if (seed_index >= num_seeds) { // just in case.
              for (int i = 0; i < num_seeds; ++i) {
                seed_values[i] = rndgen2.random_number(INT_MAX);
              }
              seed_index = 0;
            }
          }
          std::discrete_distribution pick_particle(weights.begin(), weights.end());

          for (unsigned i = r.begin(); i < r.end(); ++i) {
            found_particles[i] = current_sample[pick_particle(rndgen2.rndgen_)];
            found_particles[i].perturb(rndgen2);
            found_particles[i].sim(crown_age, min_lin, max_lin);
            if (found_particles[i].success == true) {
              double dist = calc_dist(found_particles[i]);
              if (dist < threshold[iteration]) {
                found_particles[i].update_weight(current_sample, rndgen2);
              } else {
                found_particles[i].success = false;
              }
            }
          }
      });

      int num_accepted = 0;
      for (const auto& i : found_particles) {
        if (i.success && new_sample.size() < num_particles) {
          new_sample.push_back(i);
          num_accepted++;
        }
      }
      num_tried += loop_size;
      accept_rate = 1.0 * (1 + num_accepted) / num_tried;
      Rcpp::Rcout << "\ncurrent_accept_rate: " << accept_rate << "\n";
      if (accept_rate < limit_accept_rate) {
        break;
      }
    }

    current_sample = new_sample;
    double sum_weights = 0.0;
    for (const auto& i : current_sample) {
      sum_weights += i.weight;
    }
    double factor = 1.0 / sum_weights;
    for (auto& i : current_sample) {
      i.weight *= factor;
    }

    Rcpp::Rcout << "\n";
    Rcpp::checkUserInterrupt();
  }


  double calc_dist(const particle_par& p) {
    double d1 = (p.gamma - ref_gamma);
    double d2 = (p.colless - ref_colless);
    double d3 = (p.num_lin- ref_num_lin);

    static double abs_ref_gamma = std::abs(ref_gamma);

    std::array<double, 3> diff;
    diff[0] = (d1 * d1) * 1.0 / abs_ref_gamma;
    diff[1] = (d2 * d2) * 1.0 / ref_colless;
    diff[2] = (d3 * d3) * 1.0 / ref_num_lin;

    return std::accumulate(diff.begin(), diff.end(), 0.0);   //*std::max_element(diff.begin(), diff.end());
  }

  void update_kernel(size_t iter) {
      // each sigma is to be updated to 2*var(parval)

      // calculate means
      std::array< double, 5 > means = {0.0};
      for (const auto& i : current_sample) {
        for (size_t j = 0; j < 5; ++j) {
          means[j] += std::log(i.params_[j]);
        }
      }
      for (size_t j = 0; j < 5; ++j) {
        means[j] *= 1.0 / current_sample.size();
      }

      sigmas = {0.0, 0.0, 0.0, 0.0, 0.0};
      for (const auto& i : current_sample) {
        for (size_t j = 0; j < 5; ++j) {
          double d = std::log(i.params_[j]) - means[j];
          sigmas[j] += d * d;
        }
      }
      std::cerr << "iteration:" << iter << " ";
      for (size_t j = 0; j < 5; ++j) {
        sigmas[j] *= 2.0 / current_sample.size(); // 2 VAR!
        std::cerr << sigmas[j] << " ";
      } std::cerr << "\n";
      rndgen_.update_sigmas(sigmas);
  }
};
