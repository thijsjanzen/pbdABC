#pragma once
#include <array>
#include <cmath>
#include <limits>
#include <random>

// 5x5 matrix, used for the lower-triangular Cholesky factor of the
// perturbation kernel's covariance matrix
using kernel_matrix = std::array<std::array<double, 5>, 5>;

struct rnd_t {
  std::mt19937_64 rndgen_;

  rnd_t() {
    std::random_device rd;
    rndgen_ = std::mt19937_64(rd());
    unif_dist = std::uniform_real_distribution<>(0, 1.0);
  }

  rnd_t(size_t seed) {
    rndgen_ = std::mt19937_64(seed);
    unif_dist = std::uniform_real_distribution<>(0, 1.0);
  }

  rnd_t(size_t seed,
        const rnd_t& other) {
    kernel_L = other.kernel_L;
    kernel_log_det_L = other.kernel_log_det_L;
    rndgen_ = std::mt19937_64(seed);
  }

  double uniform() {
    return unif_dist(rndgen_);
  }

  double uniform(double lower, double upper) {
    std::uniform_real_distribution<double> d(lower, upper);
    return d(rndgen_);
  }

  int random_number(unsigned int n) {
    return std::uniform_int_distribution<> (0, n-1)(rndgen_);
  }

  double exp(double lambda) {
    if (lambda == 0.0) return 1e20f;
    return std::exponential_distribution<double>(lambda)(rndgen_);
  }

  double normal(double m, double s) {
    std::normal_distribution<double> d(m, s);
    return(d(rndgen_));
  }

  // ---- multivariate normal perturbation kernel on log-parameters ----
  // The kernel covariance is C = L L^T, with L lower triangular (Cholesky).

  static kernel_matrix identity_kernel() {
    kernel_matrix I{};
    for (size_t i = 0; i < 5; ++i) I[i][i] = 1.0;
    return I;
  }

  void update_kernel(const kernel_matrix& L) {
    kernel_L = L;
    kernel_log_det_L = 0.0;
    for (size_t i = 0; i < 5; ++i) kernel_log_det_L += std::log(L[i][i]);
  }

  // correlated step L z, with z ~ N(0, I), so the step ~ N(0, L L^T)
  std::array<double, 5> kernel_step() {
    std::array<double, 5> z;
    for (auto& v : z) v = normal(0.0, 1.0);
    std::array<double, 5> step{};
    for (size_t i = 0; i < 5; ++i) {
      for (size_t j = 0; j <= i; ++j) step[i] += kernel_L[i][j] * z[j];
    }
    return step;
  }

  // log density of step d under N(0, L L^T):
  // solve L y = d (forward substitution), then -0.5 |y|^2 - log|L| - 2.5 log(2 pi)
  double kernel_log_density(const std::array<double, 5>& d) const {
    std::array<double, 5> y;
    double q = 0.0;
    for (size_t i = 0; i < 5; ++i) {
      double s = d[i];
      for (size_t j = 0; j < i; ++j) s -= kernel_L[i][j] * y[j];
      y[i] = s / kernel_L[i][i];
      q += y[i] * y[i];
    }
    return -0.5 * q - kernel_log_det_L - 2.5 * std::log(2.0 * 3.141592653589793238);
  }

  // standard deviation of the kernel along each parameter (for reporting)
  std::array<double, 5> kernel_sds() const {
    std::array<double, 5> sd{};
    for (size_t i = 0; i < 5; ++i) {
      double v = 0.0;
      for (size_t j = 0; j <= i; ++j) v += kernel_L[i][j] * kernel_L[i][j];
      sd[i] = std::sqrt(v);
    }
    return sd;
  }

  kernel_matrix kernel_L = identity_kernel();
  double kernel_log_det_L = 0.0;  // sum_i log L_ii = 0.5 * log det(C)
  std::uniform_real_distribution<> unif_dist;
};
