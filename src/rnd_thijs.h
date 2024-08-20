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
        const std::vector<double>& lambda_vals,
        double s) :
    lambdas(lambda_vals) {
    rndgen_ = std::mt19937_64(seed);
    unif_dist = std::uniform_real_distribution<>(0, 1.0);
    perturb_kernel = std::normal_distribution<double>(0.0, sigma);
    for (size_t i = 0; i < lambdas.size(); ++i) {
      log_lambdas.push_back(std::log(lambdas[i]));
    }
    sigma = s;
  }

  rnd_t(const std::vector<double>& lambda_vals,
        double s) :
    lambdas(lambda_vals) {
    std::random_device rd;
    rndgen_ = std::mt19937_64(rd());
    unif_dist = std::uniform_real_distribution<>(0, 1.0);
    perturb_kernel = std::normal_distribution<double>(0.0, sigma);
    for (size_t i = 0; i < lambdas.size(); ++i) {
      log_lambdas.push_back(std::log(lambdas[i]));
    }
    sigma = s;
  }

  std::uniform_real_distribution<> unif_dist;

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

  double perturb_particle_val(double m) {
    double new_val = std::exp( log(m) + perturb_kernel(rndgen_));
    if (isinf(new_val)) new_val = 1e12;
    return new_val;
  }

  std::array<double, 5> draw_from_prior() {
    std::array<double, 5> out;
    for (size_t i = 0; i < out.size(); ++i) {
      out[i] = exp(lambdas[i]);
    }
    return out;
  }

  double dens_prior(const std::array<double, 5>& params) const {
    double p = 0.0;
    for (size_t i = 0; i < params.size(); ++i) {
      if (params[i] < 0) return 0.0;

      p += log_lambdas[i] - lambdas[i] * params[i];
    }
    return std::exp(p);
  }

  std::vector<double> lambdas;
  std::vector<double> log_lambdas;
  double sigma;

  std::normal_distribution<double> perturb_kernel;
};
