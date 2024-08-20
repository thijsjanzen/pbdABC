//
//  mutable_dist.h
//  test_pbd2
//
//  Created by thijsjanzen on 20/08/2024.
//  // adopted from Hanno's rndutils

// notice that this does not have failsafes for all zero!

#pragma once

#include <vector>
#include <random>

struct fast_dist{
  template <typename InIt>
  fast_dist(InIt first, InIt last)
  {
    mutate(first, last);
  }

  fast_dist()
    : cdf_(1, double{ 0 })
  { // default ctor
  }

  double sum_val() {
    return cdf_.back();
  }

  template <typename Reng>
  size_t operator()(Reng& reng) const
  {
    std::uniform_real_distribution<double> d(0, 1.0);
    auto p = cdf_.back() * d(reng);
    return static_cast<size_t>(std::distance(cdf_.cbegin(), std::lower_bound(cdf_.cbegin(), cdf_.cend(), p)));
  }

  template <typename InIt>
  void mutate(InIt first, InIt last)
  {
    mutate_transform_partial(first, last, 0);
  }

  template <typename InIt>
  void mutate_transform_partial(InIt first, InIt last, size_t ofs)
  {
    const auto N = static_cast<size_t>(std::distance(first, last));
    mutate_transform_partial_n(first, N, ofs);
  }

  template <typename InIt>
  void mutate_transform_partial_n(InIt first, size_t N, size_t ofs)
  {
    if (N == 0) {
      // degenerated distribution is fine but special case
      cdf_.assign(1, double(1));
      return;
    }
    cdf_.resize(N);
    double sum(ofs > 0 ? cdf_[ofs - 1] : double(0));
    for (size_t i = ofs; i < N; ++i, ++first) {
      auto x = static_cast<double>(*first);
      cdf_[i] = sum += x;
    }
  }

  void append(double val) {
    cdf_.push_back(cdf_.back() + val);
  }

private:
  std::vector<double> cdf_;
};
