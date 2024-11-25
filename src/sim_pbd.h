#pragma once

#include "util.h"
#include "ltable.h"
#include "fast_dist.h"

#include <algorithm>
#include <vector>
#include <array>
#include <random>
#include <string>

#include <set>

enum species_status {good, incipient, extinct};
enum event {speciation, extinction, completion};

template< class ForwardIt >
size_t draw_index(ForwardIt first, ForwardIt last, std::mt19937_64& rndgen,
                  double max_val = -1) {

  if (max_val < 0.0) max_val = *std::max_element(first, last);

  auto n = std::distance(first, last);

  std::uniform_int_distribution<int> d(0, n - 1);
  std::uniform_real_distribution<double> u(0.0, 1.0);

  size_t index = 0;
  while (true) {
    index = d(rndgen);
    auto f = *(first + index) / max_val;
    if (f >= 1.0) break;
    if (u(rndgen) < f) break;
  }
  return index;
}

struct sim_pbd {
  const std::array<double, 2> spec_rate;
  const std::array<double, 2> ext_rate;
  const double compl_rate;

  std::vector< species_status > pop;
  std::vector< double > pop_spec;
  std::vector< double > pop_ext;
  std::vector< double > pop_compl;
  std::vector< size_t > pop_sp_number;

  fast_dist fast_dist_spec;
  fast_dist fast_dist_ext;
  fast_dist fast_dist_compl;

  double total_rate;
  double t;
  double max_t;
  size_t sp_number;

  int num_good_species;
  int num_incipient_species;
  std::string status;
  size_t upper_limit_species;

  std::array<size_t, 2> crown_count;

  ltable L;

  std::mt19937_64 rndgen_;

  sim_pbd(double lambda0,
          double lambda1,
          double mu0,
          double mu1,
          double com_rate,
          double maximum_time,
          int max_num) :
    spec_rate{lambda0, lambda1},
    ext_rate{mu0, mu1},
    compl_rate(com_rate),
    max_t(maximum_time),
    upper_limit_species(max_num) {
      t = 0.0;

      std::random_device d;
      std::mt19937_64 rndgen_t(d());
      rndgen_ = rndgen_t;
    }

  sim_pbd(const std::array<double, 5>& p,
          double maximum_time,
          int max_num) :
    spec_rate{p[0], p[1]},
    ext_rate{p[2], p[3]},
    compl_rate(p[4]),
    max_t(maximum_time),
    upper_limit_species(max_num) {
      t = 0.0;

      std::random_device d;
      std::mt19937_64 rndgen_t(d());
      rndgen_ = rndgen_t;
    }

  void clear_containers() {
    L.clear();
    pop.clear();
    pop_spec.clear();
    pop_ext.clear();
    pop_compl.clear();
    pop_sp_number.clear();
  }

  void run() {
    t = 0.0;
    sp_number = 1;

    clear_containers();

    L.push_back({0.0, 0, -1, -1});
    L.push_back({0.0, -1, 2, -1});

    pop.push_back(good);
    pop.push_back(good);
    pop_spec.push_back(spec_rate[species_status::good]);
    pop_spec.push_back(spec_rate[species_status::good]);
    pop_ext.push_back(ext_rate[species_status::good]);
    pop_ext.push_back(ext_rate[species_status::good]);
    pop_compl.push_back(0.0);
    pop_compl.push_back(0.0);

    pop_sp_number.push_back(sp_number);
    pop_sp_number.push_back(++sp_number);

    fast_dist_spec  = fast_dist(pop_spec.begin(), pop_spec.end());
    fast_dist_ext   = fast_dist(pop_ext.begin(), pop_ext.end());
    fast_dist_compl = fast_dist(pop_compl.begin(), pop_compl.end());

    crown_count = {1, 1};

    num_good_species = 2;
    num_incipient_species = 0;

    status = "success";
    while (true) {
      update_rates();
      double dt = draw_dt();
      t += dt;
      if (t > max_t) break;

      auto event = draw_event();
      apply_event(event);
      if (num_good_species + num_incipient_species < 2) {
        status = "extinct";
        break; // extinction
      }
      if (sp_number > upper_limit_species) {
        status = "overflow";
        break; // overflow
      }

      if (num_incipient_species > 1e5) {
        status = "overflow_incipient";
        break; // overflow
      }

      if (crown_count[0] < 1 || crown_count[1] < 1) {
        status = "extinct";
        break; // extinction
      }
    }

    // now we modify the ltable
    if (status == "success") {
      select_random();
      correct_ltable();
    }

    return;
  }

  void update_rates() {
    total_rate = fast_dist_ext.sum_val() +
      fast_dist_spec.sum_val() +
      fast_dist_compl.sum_val();
  }

  double draw_dt() {
    std::exponential_distribution<double> exp_dist(total_rate);
    return exp_dist(rndgen_);
  }

  event draw_event() {
    std::uniform_real_distribution<double> d(0.0, total_rate);
    auto r = d(rndgen_);
    r -= fast_dist_compl.sum_val();
    if (r <= 0.0) return completion;

    r -= fast_dist_spec.sum_val();
    if (r <= 0.0) return speciation;

    return extinction;
  }

  void apply_event(const event& e) {
    if (e == completion) {
      do_completion();
    }
    if (e == speciation) {
      do_speciation();
    }
    if (e == extinction) {
      do_extinction();
    }
  }

  void do_extinction() {
    auto index = fast_dist_ext(rndgen_);

    if (pop[index] == species_status::good) num_good_species--;
    if (pop[index] == species_status::incipient) num_incipient_species--;

    pop_ext[index] = 0.0;
    pop_spec[index] = 0.0;
    pop_compl[index] = 0.0;

    pop[index] = extinct;

    L[index][species_property::death_time] = t;

    if (L[index][species_property::id] < 0) {
      crown_count[0]--;
    } else {
      crown_count[1]--;
    }

    fast_dist_ext.mutate_transform_partial(pop_ext.begin() + index, pop_ext.end(), index);
    fast_dist_spec.mutate_transform_partial(pop_spec.begin() + index, pop_spec.end(), index);
    fast_dist_compl.mutate_transform_partial(pop_compl.begin() + index, pop_compl.end(), index);
  }

  void do_completion() {
    size_t index = fast_dist_compl(rndgen_);


    pop_spec[index] = spec_rate[species_status::good];
    pop_ext[index] = ext_rate[species_status::good];
    pop_compl[index] = 0.0;
    pop[index] = species_status::good;
    pop_sp_number[index] = sp_number;
    sp_number++;
    num_good_species++;
    num_incipient_species--;


    fast_dist_ext.mutate_transform_partial_n(pop_ext.begin() + index, pop_ext.size(), index);
    fast_dist_spec.mutate_transform_partial_n(pop_spec.begin() + index, pop_spec.size(), index);
    fast_dist_compl.mutate_transform_partial_n(pop_compl.begin() + index, pop_compl.size(), index);
  }

  void do_speciation() {
    //
    size_t index =  fast_dist_spec(rndgen_);


    pop_spec.push_back(spec_rate[species_status::incipient]);
    pop_ext.push_back(ext_rate[species_status::incipient]);
    pop_compl.push_back(compl_rate);
    pop.push_back(species_status::incipient);
    pop_sp_number.push_back(pop_sp_number[index]);

    // amend L table
    auto p_id = L[index][species_property::id];
    auto new_id = 1 + std::abs(L.back()[species_property::id]);
    if (p_id < 0) new_id *= -1;
    L.push_back({t, p_id, new_id, -1});
    num_incipient_species++;

    if (L.back()[species_property::id] < 0) {
      crown_count[0]++;
    } else {
      crown_count[1]++;
    }

    fast_dist_ext.append(pop_ext.back());
    fast_dist_spec.append(pop_spec.back());
    fast_dist_compl.append(pop_compl.back());
  }

  void select_random() {
    std::vector<size_t> indices(L.size());
    std::iota(indices.begin(), indices.end(), 0);
    std::shuffle(indices.begin(), indices.end(), rndgen_);

    std::set< size_t > species;

    double artificial_death_time = max_t - 1e-3;

    for (const auto& i : indices) {
      if (species.find(pop_sp_number[i]) != species.end()) {
        L[i][species_property::death_time] = artificial_death_time;
      } else {
        species.insert(pop_sp_number[i]);
      }
    }
    return;
  }

  void correct_ltable() {
    // make it go from crown age to 0 instead of the other way around
    for (auto& i : L) {
      i[species_property::birth_time] = max_t - i[species_property::birth_time];

      if (i[species_property::death_time] >= 0.0) {
        i[species_property::death_time] = max_t - i[species_property::death_time];
      }
    }
  }

  size_t get_num_good_species() {
    return num_good_species;
  }

  size_t get_num_incipient_species() {
    return num_incipient_species;
  }
};

inline std::vector< std::array<double, 4>> sim_once(const std::array<double, 5>& params,
                                                    double crown_age,
                                                    int max_num_species,
                                                    bool *success,
                                                    int *num_lin) {

  sim_pbd sim(params, crown_age, max_num_species);
  *success = false;
  sim.run();

  if (sim.status == "success") *success = true;

  *num_lin = 0;
  for (const auto& i : sim.L) {
    if (i[species_property::death_time] == -1) (*num_lin)++;
  }

  return(sim.L);
}

