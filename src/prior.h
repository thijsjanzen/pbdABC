#pragma once

#include <cmath>
#include "rnd_thijs.h"

enum dist_type {expon, uniform, inv_expon, nothing};

const size_t num_params = 5;

using param_set = std::array<double, num_params>;

struct prior {

  prior() {
    prior_type = dist_type::nothing;
  }

  prior(const param_set& lambdas) :
    means(lambdas) {
    prior_type = dist_type::expon;
    for (size_t i = 0; i < lambdas.size(); ++i) {
      log_means[i] = std::log(lambdas[i]);
    }
  }

  prior(const param_set& lambdas,
        bool inv_compl_rate) :
    means(lambdas) {
    prior_type = dist_type::expon;
    for (size_t i = 0; i < lambdas.size(); ++i) {
      log_means[i] = std::log(lambdas[i]);
    }
    if (inv_compl_rate == true) {
   //    log_means.back() = std::log(1.0 / lambdas.back());
   //    means.back() = 1.0 / lambdas.back();
       prior_type = dist_type::inv_expon;
    }
  }

  prior(const param_set& low,
        const param_set& up) :
    lower(low),
    upper(up) {
    prior_type = dist_type::uniform;
  }

  param_set lower;
  param_set upper;
  param_set means;
  param_set log_means;

  dist_type prior_type;

  double dens_prior(const param_set& p) const {
    //return prior_type == dist_type::uniform ? dens_uniform(p) : dens_exp(p);
    if (prior_type == dist_type::uniform) return dens_uniform(p);
    if (prior_type == dist_type::expon) return dens_exp(p);
    if (prior_type == dist_type::inv_expon) return dens_inv_exp(p);

    return dens_exp(p); // default
  }

  bool pass_prior(const param_set& p) const {
    if (prior_type == dist_type::uniform) return pass_uniform(p);
    if (prior_type == dist_type::expon) return pass_exp(p);
    if (prior_type == dist_type::inv_expon) return pass_inv_exp(p);

    // return prior_type == dist_type::uniform ? pass_uniform(p) : pass_exp(p);
    //
    return pass_exp(p); // default
  }

  param_set gen_prior(rnd_t& rndgen) const {
    if (prior_type == dist_type::uniform) return gen_uniform(rndgen);
    if (prior_type == dist_type::expon) return gen_exp(rndgen);
    if (prior_type == dist_type::inv_expon) return gen_inv_exp(rndgen);

   // return prior_type == dist_type::uniform ? gen_uniform(rndgen) : gen_exp(rndgen);
    return gen_exp(rndgen); // default
  }

  double dens_uniform(const param_set& p) const {
    for (size_t i = 0; i < num_params; ++i) {
      if (std::log10(p[i]) < lower[i]) return 0.0;
      if (std::log10(p[i]) > upper[i]) return 0.0;
    }
    return 1.0;
  }

  double dens_exp(const param_set& params) const {
    double p = 0.0;
    for (size_t i = 0; i < num_params; ++i) {
      if (params[i] < 0) return 0.0;

      p += log_means[i] - means[i] * params[i];
    }
    return std::exp(p);
  }

  double dens_inv_exp(const param_set& params) const {
    double p = 0.0;
    for (size_t i = 0; i < num_params - 1; ++i) {
      if (params[i] < 0) return 0.0;

      p += log_means[i] - means[i] * params[i];
    }

    if (params.back() < 0) return 0.0;
    p += log_means.back() - means.back() * (1.0 / params.back());

    return std::exp(p);
  }

  double dens_inv_only_compl_exp(const param_set& params) const {
    double p = log_means.back() - means.back() * (1.0 / params.back());
    return std::exp(p);
  }

  double dens_only_compl_exp(const param_set& params) const {
    double p = log_means.back() - means.back() * (params.back());
    return std::exp(p);
  }


  bool pass_uniform(const param_set& params) const {
    auto prob_prior = dens_prior(params);
    if (prob_prior > 0.0) return true;

    return false;
  }

  bool pass_exp(const param_set& params) const {
    for (const auto& i : params) {
      if (i < 0.0) return false;
    }
    return true;
  }

  bool pass_inv_exp(const param_set& params) const {
    for (const auto& i : params) {
      if (i < 0.0) return false;
    }
    return true;
  }

  param_set gen_uniform(rnd_t& rndgen) const {
    param_set out;
    for (size_t i = 0; i < num_params; ++i) {
      out[i] = pow(10, rndgen.uniform(lower[i], upper[i]));
    }
    return out;
  }

  param_set gen_exp(rnd_t& rndgen) const {
    param_set out;
    for (size_t i = 0; i < num_params; ++i) {
      out[i] = rndgen.exp(means[i]);
    }
    return out;
  }

  param_set gen_inv_exp(rnd_t& rndgen) const {
    param_set out;
    for (size_t i = 0; i < num_params; ++i) {
      out[i] = rndgen.exp(means[i]);
    }
    out.back() = 1.0 / out.back();

    return out;
  }
};
