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
#include <chrono>
#include <cstdio>
#include <string>
#include <cstdint>
#include <random>

#include <RcppParallel.h>

// Non-throwing check for a user interrupt (stop button / Esc / Ctrl-C).
// Rcpp::checkUserInterrupt() throws, which unwinds everything and loses the
// results; R_ToplevelExec runs the check in its own top-level context and
// returns FALSE if it was interrupted. Must be called from the main thread.
// Colless index corrected for tree size under the Yule model, as
// treestats::colless(phy, normalization = "yule"): (C - n log n -
// n (gamma_E - 1 - log 2)) / n. Unlike the raw index it hardly depends on n.
inline double colless_yule(double colless, double num_tips) {
  const double g = 0.577215664901532;  // Euler-Mascheroni constant
  return (colless - num_tips * std::log(num_tips) -
          num_tips * (g - 1.0 - std::log(2.0))) / num_tips;
}

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
  const int num_generations;
  const double limit_accept_rate;
  // if in (0, 1): adaptive threshold, i.e. the acceptance threshold of each
  // generation is this quantile of the distances of the previous generation.
  // Otherwise the fixed schedule 10 * exp(-0.5 * (i - 1)) is used.
  const double threshold_quantile;
  // perturbation kernel covariance = kernel_scale * weighted covariance of
  // the log-parameters of the current population
  const double kernel_scale;
  // scaled distance: sum over statistics of ((sim - obs) / scale)^2, with
  // the scale of each statistic its spread (MAD) in generation 0 and the
  // Colless index normalised for tree size. Otherwise the original relative
  // distance sum (sim - obs)^2 / |obs| with the raw Colless index.
  const bool scaled_distance;
  std::array<double, 3> stat_scale = {1.0, 1.0, 1.0};
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
           double k_scale = 2.0,
           bool scaled_dist = false) :
    ref_gamma(obs_gamma),
    ref_colless(obs_colless),
    ref_num_lin(obs_num_lin),
    crown_age(ca),
    min_lin(minimum_lineages),
    max_lin(maximum_lineages),
    num_particles(n),
    num_generations(num_iterations),
    limit_accept_rate(limit_rate),
    threshold_quantile(thresh_quantile),
    kernel_scale(k_scale),
    scaled_distance(scaled_dist) {
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

  // One progress line per generation, redrawn in place after every batch of
  // simulations and kept (newline) once the generation is finished.
  size_t last_progress_length = 0;
  std::chrono::steady_clock::time_point generation_start;

  void print_progress(int generation, size_t done, long num_found,
                      long num_tried, bool show_threshold, bool finished) {
    const int width = 30;
    const size_t total = static_cast<size_t>(num_particles);
    const size_t shown = std::min(done, total);
    const int filled = total > 0 ? static_cast<int>(width * shown / total)
                                 : width;
    std::string bar(filled, '=');
    if (filled < width) {
      bar += '>';
      bar += std::string(width - filled - 1, ' ');
    }
    const int count_width = static_cast<int>(std::to_string(total).size());
    const double seconds = std::chrono::duration<double>(
      std::chrono::steady_clock::now() - generation_start).count();

    std::string line;
    char buf[160];
    std::snprintf(buf, sizeof(buf), "Generation %*d/%d [%s] %*zu/%zu",
                  static_cast<int>(std::to_string(num_generations - 1).size()),
                  generation, num_generations - 1, bar.c_str(),
                  count_width, shown, total);
    line += buf;
    if (num_tried > 0) {
      std::snprintf(buf, sizeof(buf), " | acceptance %.3g%%",
                    100.0 * num_found / num_tried);
      line += buf;
    }
    if (show_threshold) {
      std::snprintf(buf, sizeof(buf), " | threshold %.4g", current_threshold);
      line += buf;
    }
    std::snprintf(buf, sizeof(buf), " | %.0f s", seconds);
    line += buf;

    // overwrite leftovers of a longer previous line
    const size_t length = line.size();
    if (length < last_progress_length) {
      line += std::string(last_progress_length - length, ' ');
    }
    last_progress_length = finished ? 0 : length;

    Rcpp::Rcout << "\r" << line;
    if (finished) Rcpp::Rcout << "\n";
    Rcpp::Rcout << std::flush;
  }

  void iterate_first() {
    accept_rate = 1.0;
    long num_tried = 0;
    long num_found = 0;
    generation_start = std::chrono::steady_clock::now();
    print_progress(0, current_sample.size(), num_found, num_tried, false, false);

    while(current_sample.size() < num_particles) {
      if (user_interrupted()) {
        interrupted = true;
        return;
      }

      const int num_remaining = num_particles - current_sample.size();
      int loop_size = num_remaining * 1.0 / accept_rate;
      if (loop_size > 10000) loop_size = 10000;
      if (loop_size < 0) loop_size = 10000;

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

      for (const auto& i : found_particles) {
        if (i.success) {
          num_found++;
          if (current_sample.size() < num_particles) current_sample.push_back(i);
        }
      }
      num_tried += loop_size;
      accept_rate = 1.0 * (1 + num_found) / num_tried;
      print_progress(0, current_sample.size(), num_found, num_tried, false,
                     false);
    }

    print_progress(0, current_sample.size(), num_found, num_tried, false, true);
    if (scaled_distance) compute_stat_scales();
    update_distances();
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
    long num_tried = 0;
    long num_found = 0;  // particles accepted this generation (incl. surplus)
    generation_start = std::chrono::steady_clock::now();
    print_progress(iteration, 0, num_found, num_tried, true, false);

    while(new_sample.size() < num_particles) {
      if (user_interrupted()) {
        // discard the unfinished generation; current_sample is unchanged
        interrupted = true;
        return;
      }

      num_remaining = num_particles - new_sample.size();

      int loop_size = num_remaining * 1.0 / accept_rate;
      // batches are capped so that a stop is noticed within one batch
      if (loop_size > 10000) loop_size = 10000;
      if (loop_size < 0) loop_size = 10000;

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
      print_progress(iteration, new_sample.size(), num_found, num_tried, true,
                     false);
      if (accept_rate < limit_accept_rate) {
        break;
      }
    }
    print_progress(iteration, new_sample.size(), num_found, num_tried, true,
                   true);
    if (new_sample.size() < num_particles) {
      Rcpp::Rcout << "acceptance rate below limiting_accept_rate, stopping\n";
    }

    if (new_sample.empty()) {
      // nothing accepted: keep the previous population, caller stops
      Rcpp::Rcout << "no particles accepted\n";
      stopped_early = true;
      return;
    }

    current_sample = new_sample;
    double sum_weights = 0.0;
    for (const auto& i : current_sample) {
      sum_weights += i.weight;
    }
    if (!(sum_weights > 0.0) || !std::isfinite(sum_weights)) {
      Rcpp::Rcout << "invalid particle weights, stopping\n";
      stopped_early = true;
      return;
    }

    double factor = 1.0 / sum_weights;
    for (auto& i : current_sample) {
      i.weight *= factor;
    }
  }


  // the summary statistics as used by the scaled distance
  std::array<double, 3> distance_stats(double gamma, double colless,
                                       double num_lin) const {
    return {gamma, colless_yule(colless, num_lin), num_lin};
  }

  // scale of each statistic: MAD (x 1.4826, a robust standard deviation) in
  // the current population, i.e. generation 0; falls back to the standard
  // deviation, and to 1, if a statistic does not vary
  void compute_stat_scales() {
    const size_t n = current_sample.size();
    for (size_t k = 0; k < 3; ++k) {
      std::vector<double> x;
      x.reserve(n);
      for (const auto& p : current_sample) {
        x.push_back(distance_stats(p.gamma, p.colless, p.num_lin)[k]);
      }
      if (x.size() < 2) { stat_scale[k] = 1.0; continue; }
      auto median_of = [](std::vector<double> v) {
        std::nth_element(v.begin(), v.begin() + v.size() / 2, v.end());
        return v[v.size() / 2];
      };
      const double med = median_of(x);
      std::vector<double> dev(x.size());
      for (size_t i = 0; i < x.size(); ++i) dev[i] = std::abs(x[i] - med);
      double scale = 1.4826 * median_of(dev);
      if (!(scale > 0.0) || !std::isfinite(scale)) {
        double mean = 0.0, var = 0.0;
        for (double v : x) mean += v / x.size();
        for (double v : x) var += (v - mean) * (v - mean) / (x.size() - 1);
        scale = std::sqrt(var);
      }
      stat_scale[k] = (scale > 0.0 && std::isfinite(scale)) ? scale : 1.0;
    }
    Rcpp::Rcout << "distance scales (MAD in generation 0): gamma "
                << stat_scale[0] << ", Colless (size-normalised) "
                << stat_scale[1] << ", number of species " << stat_scale[2]
                << "\n";
  }

  // the scales define the distance, so a resumed run must reuse them
  void write_scales(int sim_number) const {
    std::ofstream out("scales_" + std::to_string(sim_number) + ".txt");
    out.precision(17);
    out << stat_scale[0] << " " << stat_scale[1] << " " << stat_scale[2] << "\n";
  }

  bool read_scales(int sim_number) {
    std::ifstream in("scales_" + std::to_string(sim_number) + ".txt");
    std::array<double, 3> s;
    if (!(in >> s[0] >> s[1] >> s[2])) return false;
    stat_scale = s;
    return true;
  }

  double calc_dist(const particle& p) {
    if (scaled_distance) {
      const auto sim = distance_stats(p.gamma, p.colless, p.num_lin);
      const auto obs = distance_stats(ref_gamma, ref_colless, ref_num_lin);
      double d = 0.0;
      for (size_t k = 0; k < 3; ++k) {
        const double z = (sim[k] - obs[k]) / stat_scale[k];
        d += z * z;
      }
      return std::isfinite(d) ? d : std::numeric_limits<double>::infinity();
    }

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

        Rcpp::Rcout << "found previous output up to generation " << iter
                    << ", continuing from there\n";

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
        if (scaled_distance && !read_scales(sim_number)) {
          // output from a run without a scales file: estimate them from the
          // resumed population (the distance then differs from that run)
          Rcpp::Rcout << "no scales_" << sim_number << ".txt found; "
                      << "estimating distance scales from the resumed "
                      << "population\n";
          compute_stat_scales();
        }
        update_distances();
        starting_iteration = iter + 1;
        blank_sheet = false;

        Rcpp::Rcout << "read " << current_sample.size() << " particles\n";

        break;
      }
    }
  }
};
