#pragma once

#include "prior.h"
#include "rnd_thijs.h"
#include "sim_pbd.h"
#include "ltable.h"

#include "L2newick.h"

#include <limits>

struct particle {
  param_set params_;
  double gamma;
  double colless;
  int num_lin;

  bool success = false;

  double weight = 1.0;
  double dist = std::numeric_limits<double>::max();  // distance to the data
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

  // Perturb all parameters jointly on the log scale with a correlated
  // multivariate normal step (covariance = kernel_L kernel_L^T).
  void perturb(rnd_t& rndgen) {
    auto step = rndgen.kernel_step();
    for (size_t i = 0; i < params_.size(); ++i) {
      double new_val = std::exp(std::log(params_[i]) + step[i]);
      if (std::isinf(new_val)) new_val = 1e12;
      // zero is deadly, let's use small, but not toooo small.
      if (new_val <= 0.0) new_val = 10 * std::numeric_limits<double>::min();
      params_[i] = new_val;
    }
    return;
  }

  double prob_perturb(const particle& other,
                      const rnd_t& rndgen) {

    // Density of moving from 'other' to this particle in log-parameter space,
    // under the multivariate normal kernel.
    std::array<double, 5> d;
    for (size_t i = 0; i < 5; ++i) {
      d[i] = std::log(params_[i]) - std::log(other.params_[i]);
    }
    return std::exp(rndgen.kernel_log_density(d));
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

// Cholesky factorisation of a symmetric 5x5 matrix A + jitter * I.
// Returns false if the matrix is not (numerically) positive definite.
inline bool cholesky_5(const kernel_matrix& A, double jitter, kernel_matrix* L) {
  kernel_matrix out{};
  for (size_t i = 0; i < 5; ++i) {
    for (size_t j = 0; j <= i; ++j) {
      double s = A[i][j] + (i == j ? jitter : 0.0);
      for (size_t k = 0; k < j; ++k) s -= out[i][k] * out[j][k];
      if (i == j) {
        if (!(s > 0.0) || !std::isfinite(s)) return false;
        out[i][i] = std::sqrt(s);
      } else {
        out[i][j] = s / out[j][j];
      }
    }
  }
  *L = out;
  return true;
}

// Multivariate normal kernel with covariance scale * (weighted covariance of
// the log-parameters of the current population); scale = 2 follows Beaumont
// et al. (2009) / Filippi et al. (2013). Returns the Cholesky factor.
inline kernel_matrix compute_kernel(const std::vector<particle>& sample,
                                    double scale = 2.0) {
  kernel_matrix L = rnd_t::identity_kernel();
  if (sample.size() < 2) return L;

  double sum_w = 0.0;
  for (const auto& p : sample) sum_w += p.weight;
  const bool use_weights = sum_w > 0.0 && std::isfinite(sum_w);
  auto w = [&](const particle& p) {
    return use_weights ? p.weight / sum_w : 1.0 / sample.size();
  };

  std::array<double, 5> mean{};
  for (const auto& p : sample) {
    for (size_t j = 0; j < 5; ++j) mean[j] += w(p) * std::log(p.params_[j]);
  }

  kernel_matrix cov{};
  for (const auto& p : sample) {
    std::array<double, 5> d;
    for (size_t j = 0; j < 5; ++j) d[j] = std::log(p.params_[j]) - mean[j];
    for (size_t j = 0; j < 5; ++j) {
      for (size_t k = 0; k <= j; ++k) cov[j][k] += w(p) * d[j] * d[k];
    }
  }
  double max_diag = 1e-12;
  for (size_t j = 0; j < 5; ++j) {
    for (size_t k = 0; k <= j; ++k) {
      cov[j][k] *= scale;
      cov[k][j] = cov[j][k];
      if (!std::isfinite(cov[j][k])) return L;  // fall back to identity
    }
    if (cov[j][j] < 1e-12) cov[j][j] = 1e-12;   // avoid a degenerate kernel
    max_diag = std::max(max_diag, cov[j][j]);
  }

  // a (near-)singular covariance, e.g. with few distinct particles, gets a
  // little jitter on the diagonal
  for (double jitter : {0.0, 1e-10, 1e-8, 1e-6, 1e-4, 1e-2}) {
    if (cholesky_5(cov, jitter * max_diag, &L)) return L;
  }

  // last resort: diagonal kernel with the same marginal variances
  L = kernel_matrix{};
  for (size_t j = 0; j < 5; ++j) L[j][j] = std::sqrt(cov[j][j]);
  return L;
}
