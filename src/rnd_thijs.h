#pragma once
#include <random>

struct rnd_t {
  std::mt19937_64 rndgen_;

  rnd_t() {
    std::random_device rd;
    rndgen_ = std::mt19937_64(rd());
    unif_dist = std::uniform_real_distribution<>(0, 1.0);
  }

  rnd_t(size_t seed,
        const std::vector<double>& low,
        const std::vector<double>& up) :
    lower(low),
    upper(up) {
    rndgen_ = std::mt19937_64(seed);
    unif_dist = std::uniform_real_distribution<>(0, 1.0);
  }

  rnd_t(const std::vector<double>& low,
        const std::vector<double>& up) :
    lower(low),
    upper(up) {
    std::random_device rd;
    rndgen_ = std::mt19937_64(rd());
    unif_dist = std::uniform_real_distribution<>(0, 1.0);
  }

  rnd_t(size_t seed,
        const rnd_t& other) {
    lower = other.lower;
    upper = other.upper;
    kernel_sigmas = other.kernel_sigmas;
    rndgen_ = std::mt19937_64(seed);
  }

  double uniform() {
    return unif_dist(rndgen_);
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

  double perturb_particle_val(double m, size_t i) {
    double new_val = std::exp( log(m) + normal(0.0, kernel_sigmas[i]));
    if (std::isinf(new_val)) new_val = 1e12;
    return new_val;
  }

  std::array<double, 5> draw_from_prior() {
    std::array<double, 5> out;
    for (size_t i = 0; i < out.size(); ++i) {
      std::uniform_real_distribution<double> d(lower[i], upper[i]);
      out[i] = pow(10, d(rndgen_));
    }
    return out;
  }

  double dens_prior(const std::array<double, 5>& params) const {
    for (size_t i = 0; i < params.size(); ++i) {
      if (std::log10(params[i]) < lower[i]) return 0.0;
      if (std::log10(params[i]) > upper[i]) return 0.0;
    }
    return 1.0;
  }

  void update_sigmas(std::array<double, 5> s) {
    kernel_sigmas = s;
  }

  std::vector<double> lower;
  std::vector<double> upper;

  std::array<double, 5> kernel_sigmas;
  std::uniform_real_distribution<> unif_dist;
};
