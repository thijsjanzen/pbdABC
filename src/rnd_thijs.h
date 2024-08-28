#pragma once
#include <random>

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
    kernel_sigmas = other.kernel_sigmas;
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

  double perturb_particle_val(double m, size_t i) {
    double new_val = std::exp( log(m) + normal(0.0, kernel_sigmas[i]));
    if (std::isinf(new_val)) new_val = 1e12;
    return new_val;
  }

  void update_sigmas(std::array<double, 5> s) {
    kernel_sigmas = s;
  }

  std::array<double, 5> kernel_sigmas;
  std::uniform_real_distribution<> unif_dist;
};
