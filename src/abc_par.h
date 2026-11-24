#pragma once

#include <vector>
#include <array>
#include <fstream>
#include "colless.h"
#include "gamma.h"

#include "rnd_thijs.h"
#include "sim_pbd.h"
#include "ltable.h"
#include "prior.h"
#include "particle.h"

#include <mutex>
#include <cstdint>
#include <random>

#include <RcppParallel.h>

// Non-throwing check for a user interrupt (stop button / Esc / Ctrl-C).
// Rcpp::checkUserInterrupt() throws, which unwinds everything and loses the
// results; R_ToplevelExec runs the check in its own top-level context and
// returns FALSE if it was interrupted. Must be called from the main thread.
inline void check_interrupt_fn(void* /*dummy*/) {
  R_CheckUserInterrupt();
}

inline bool user_interrupted() {
  return R_ToplevelExec(check_interrupt_fn, nullptr) == FALSE;
}


struct analysis_par {
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
  const double limit_accept_rate;
  // if in (0, 1): adaptive threshold, i.e. the acceptance threshold of each
  // generation is this quantile of the distances of the previous generation.
  // Otherwise the fixed schedule 10 * exp(-0.5 * (i - 1)) is used.
  const double threshold_quantile;
  // perturbation kernel covariance = kernel_scale * weighted covariance of
  // the log-parameters of the current population
  const double kernel_scale;
  double current_threshold = 0.0;

  double accept_rate;

  bool blank_sheet;
  bool stopped_early = false;
  bool interrupted = false;  // user pressed stop: return completed generations
  size_t starting_iteration;

  std::vector<double> threshold;

  prior prior_dist;

  analysis_par(int n,
           int num_iterations,
           double ca,
           double minimum_lineages,
           double maximum_lineages,
           param_set lower,
           param_set upper,
           param_set means,
           bool use_inv_prior,
           double obs_gamma,
           double obs_colless,
           double obs_num_lin,
           double limit_rate,
           double thresh_quantile = -1.0,
           double k_scale = 2.0) :
    ref_gamma(obs_gamma),
    ref_colless(obs_colless),
    ref_num_lin(obs_num_lin),
    crown_age(ca),
    min_lin(minimum_lineages),
    max_lin(maximum_lineages),
    num_particles(n),
    limit_accept_rate(limit_rate),
    threshold_quantile(thresh_quantile),
    kernel_scale(k_scale) {
    rndgen_ = rnd_t();
    for (size_t i = 0; i < num_iterations; ++i) {
      threshold.push_back(10 * std::exp(-0.5 * (i - 1)));
    }

    if (means.empty()) {
      prior_dist = prior(lower, upper);
    } else if (means[0] < 0) {
      prior_dist = prior(lower, upper);
    } else {
      prior_dist = prior(means, use_inv_prior);
    }
  }

  bool adaptive_threshold() const {
    return threshold_quantile > 0.0 && threshold_quantile < 1.0;
  }

  // distance of every particle in the current population to the data
  void update_distances() {
    for (auto& p : current_sample) p.dist = calc_dist(p);
  }

  // q-quantile (type 1, inverse empirical cdf) of the current distances
  double distance_quantile(double q) const {
    std::vector<double> d;
    d.reserve(current_sample.size());
    // only finite distances: NaN makes nth_element's result undefined
    for (const auto& p : current_sample) {
      if (std::isfinite(p.dist)) d.push_back(p.dist);
    }
    if (d.empty()) return std::numeric_limits<double>::max();
    size_t k = static_cast<size_t>(std::ceil(q * d.size()));
    if (k > 0) k--;
    std::nth_element(d.begin(), d.begin() + k, d.end());
    return d[k];
  }

  double threshold_for(int iteration) const {
    if (adaptive_threshold()) return distance_quantile(threshold_quantile);
    return threshold[iteration];
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


    Rcpp::Rcout << "0--------25--------50--------75--------100\n";
    Rcpp::Rcout << "*";
    int updateFreq = num_particles / 20;
    if(updateFreq < 1) updateFreq = 1;

    int prev_update = 0;

    while(current_sample.size() < num_particles) {
      if (user_interrupted()) {
        interrupted = true;
        return;
      }
      for (size_t t = prev_update; t < current_sample.size(); ++t) {
        if (t % updateFreq == 0) {
          Rcpp::Rcout << "**";
        }
      }
      prev_update = current_sample.size();


      int loop_size = num_remaining * 1.0 / accept_rate;
      if (loop_size > 10000) loop_size = 10000;
      if (loop_size < 0) loop_size = 10000;

      std::cerr << "\nloop_size: " << loop_size << " " << current_sample.size() << "\n";

      const uint64_t batch_seed = rndgen_.rndgen_();

      std::vector< particle > found_particles(loop_size);

      tbb::parallel_for(
        tbb::blocked_range<unsigned>(0, loop_size),
        [&](const tbb::blocked_range<unsigned>& r) {

          // independent stream per block: seed depends on the batch seed and
          // the block start, so no shared state between threads.
          std::seed_seq sseq{batch_seed, static_cast<uint64_t>(r.begin())};
          std::array<uint64_t, 1> block_seed;
          sseq.generate(block_seed.begin(), block_seed.end());
          rnd_t rndgen2(block_seed[0], rndgen_);

          for (unsigned i = r.begin(); i < r.end(); ++i) {

            auto new_particle = particle(prior_dist, rndgen2);

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

    update_distances();
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
    current_threshold = threshold_for(iteration);
    Rcpp::Rcout << "iteration: " << iteration
                << ", threshold: " << current_threshold << "\n";
    Rcpp::Rcout << "0--------25--------50--------75--------100\n";
    Rcpp::Rcout << "*";
    int updateFreq = num_particles / 20;
    if(updateFreq < 1) updateFreq = 1;

    int prev_update = 0;
    int num_tried = 0;
    int num_found = 0;  // particles accepted this generation (incl. surplus)
    while(new_sample.size() < num_particles) {
      if (user_interrupted()) {
        // discard the unfinished generation; current_sample is unchanged
        interrupted = true;
        return;
      }
      for (size_t t = prev_update; t < new_sample.size(); ++t) {
        if (t % updateFreq == 0) {
          Rcpp::Rcout << "**";
        }
      }
      prev_update = new_sample.size();

      num_remaining = num_particles - new_sample.size();

      int loop_size = num_remaining * 1.0 / accept_rate;
      // batches are capped so that a stop is noticed within one batch
      if (loop_size > 10000) loop_size = 10000;
      if (loop_size < 0) loop_size = 10000;

      std::cerr << "\nloop_size: " << loop_size << "\n";


      const uint64_t batch_seed = rndgen_.rndgen_();

      std::vector< particle > found_particles(loop_size);

      tbb::parallel_for(
        tbb::blocked_range<unsigned>(0, loop_size),
        [&](const tbb::blocked_range<unsigned>& r) {

          // independent stream per block: seed depends on the batch seed and
          // the block start, so no shared state between threads.
          std::seed_seq sseq{batch_seed, static_cast<uint64_t>(r.begin())};
          std::array<uint64_t, 1> block_seed;
          sseq.generate(block_seed.begin(), block_seed.end());
          rnd_t rndgen2(block_seed[0], rndgen_);
          std::discrete_distribution pick_particle(weights.begin(), weights.end());

          for (unsigned i = r.begin(); i < r.end(); ++i) {
            found_particles[i] = current_sample[pick_particle(rndgen2.rndgen_)];
            found_particles[i].perturb(rndgen2);
            if (found_particles[i].pass_prior(prior_dist)) {
              found_particles[i].sim(crown_age, min_lin, max_lin);
              if (found_particles[i].success == true) {
                double dist = calc_dist(found_particles[i]);
                found_particles[i].dist = dist;
                if (dist <= current_threshold) {
                  found_particles[i].update_weight(current_sample, rndgen2, prior_dist);
                } else {
                  found_particles[i].success = false;
                }
              }
            } else {
              found_particles[i].success = false;
            }
          }
      });

      for (const auto& i : found_particles) {
        if (i.success) {
          num_found++;
          if (new_sample.size() < num_particles) new_sample.push_back(i);
        }
      }
      // acceptance rate over the whole generation: both numerator and
      // denominator are cumulative (previously only the last batch's
      // acceptances were divided by all tries, which made the rate look far
      // lower than it was and stopped runs early)
      num_tried += loop_size;
      accept_rate = 1.0 * (1 + num_found) / num_tried;
      Rcpp::Rcout << "\ncurrent_accept_rate: " << accept_rate << "\n";
      if (accept_rate < limit_accept_rate) {
        break;
      }
    }

    if (new_sample.empty()) {
      // nothing accepted: keep the previous population, caller stops
      Rcpp::Rcout << "\nno particles accepted, stopping\n";
      stopped_early = true;
      return;
    }

    current_sample = new_sample;
    double sum_weights = 0.0;
    for (const auto& i : current_sample) {
      sum_weights += i.weight;
    }
    if (!(sum_weights > 0.0) || !std::isfinite(sum_weights)) {
      Rcpp::Rcout << "\ninvalid particle weights, stopping\n";
      stopped_early = true;
      return;
    }

    double factor = 1.0 / sum_weights;
    for (auto& i : current_sample) {
      i.weight *= factor;
    }

    Rcpp::Rcout << "\n";
  }


  double calc_dist(const particle& p) {
    double d1 = (p.gamma - ref_gamma);
    double d2 = (p.colless - ref_colless);
    double d3 = (p.num_lin- ref_num_lin);

    const double abs_ref_gamma = std::abs(ref_gamma);

    std::array<double, 3> diff;
    diff[0] = (d1 * d1) * 1.0 / abs_ref_gamma;
    diff[1] = (d2 * d2) * 1.0 / ref_colless;
    diff[2] = (d3 * d3) * 1.0 / ref_num_lin;

    const double d = std::accumulate(diff.begin(), diff.end(), 0.0);
    // a non-finite distance (e.g. from undefined statistics) never passes
    return std::isfinite(d) ? d : std::numeric_limits<double>::infinity();
  }

  void update_kernel(size_t iter) {
      // multivariate normal kernel on log-parameters with covariance
      // kernel_scale * weighted covariance of the current population
      rndgen_.update_kernel(compute_kernel(current_sample, kernel_scale));
      auto sds = rndgen_.kernel_sds();
      std::cerr << "iteration:" << iter << " kernel sd:";
      for (size_t j = 0; j < 5; ++j) {
        std::cerr << " " << sds[j];
      } std::cerr << "\n";
  }

  void write_to_file(int sim_number,
                     int iter,
                     double obs_gamma,
                     double obs_colless,
                     double obs_num_lin,
                     double bd_lambda,
                     double bd_mu) {
    std::string file_name = "res_" + std::to_string(sim_number) +
      "_" + std::to_string(iter) + ".txt";
    std::ofstream out(file_name.c_str());
    for (const auto& k : current_sample) {
      out << iter << " ";
      for (size_t j = 0; j < k.params_.size(); ++j) {
        out << k.params_[j] << " ";
      }

      out << k.gamma << " " << k.colless << " " <<
        static_cast<double>(k.num_lin) << " " << k.weight << " " <<
          obs_gamma << " " << obs_colless << " " << obs_num_lin << " "  <<
            bd_lambda << " " << bd_mu << "\n";
    }
    out.close();
    return;
  }

  void write_trees(int sim_number,
                   int iter) {
    std::string file_name = "trees_" + std::to_string(sim_number) +
      "_" + std::to_string(iter) + ".txt";
    std::ofstream out(file_name.c_str());
    for (const auto& k : current_sample) {
      out << k.newick_tree << "\n";
    }
    out.close();
  }

  void find_output(int sim_number) {
    blank_sheet = true;

    for (int iter = 100; iter >= 0; --iter) {
      std::string file_name = "res_" + std::to_string(sim_number) +
        "_" + std::to_string(iter) + ".txt";
      std::ifstream in_file(file_name.c_str());

      if (in_file.good()) {

        std::cerr << "found previous output from iteration: " << iter << "\n";

        current_sample.clear();

        particle row_entry;
        int focal_iter;
        while(!in_file.eof()) {
          in_file >>  focal_iter;
          for (size_t i = 0; i < row_entry.params_.size(); ++i) {
              in_file >> row_entry.params_[i];
          }
          // out << k.gamma << " " << k.colless << " " <<
          // static_cast<double>(k.num_lin) << " " << k.weight << " " <<
          //  obs_gamma << " " << obs_colless << " " << obs_num_lin << " "  <<
          //    bd_lambda << " " << bd_mu << "\n";
          in_file >> row_entry.gamma;
          in_file >> row_entry.colless;
          in_file >> row_entry.num_lin;
          in_file >> row_entry.weight;
          current_sample.push_back(row_entry);

          // burn other entries
          double temp;
          for (int b = 0; b < 5; ++b) {
            in_file >> temp;
          }
        }
        current_sample.pop_back(); // last entry always duplicates
        update_distances();
        starting_iteration = iter + 1;
        blank_sheet = false;

        std::cerr << "read: " << current_sample.size() << "entries\n";

        break;
      }
    }
  }
};
