//
//  fast_dist.h
//
//  Discrete distribution over lineages with mutable weights, used to pick the
//  lineage an event happens to (probability proportional to its rate).
//
//  Implemented as a Fenwick tree (binary indexed tree): slot i (1-based)
//  holds the sum of the weights in the block (i - lowbit(i), i], so that
//  changing one weight, appending a weight, the total and drawing an index
//  are all O(log n). (The previous cumulative-sum array needed O(n) per
//  changed weight, which made large simulations O(n^2).)
//

#pragma once

#include <algorithm>
#include <vector>
#include <random>
#include <cstddef>

struct fast_dist {
  fast_dist() : tree_(1, 0.0) {}

  template <typename InIt>
  fast_dist(InIt first, InIt last) {
    vals_.assign(first, last);
    rebuild();
  }

  size_t size() const { return vals_.size(); }

  // total weight
  double sum_val() const { return prefix_sum(vals_.size()); }

  // weight of index i (0-based)
  double value(size_t i) const { return vals_[i]; }

  // set the weight of index i (0-based)
  void update(size_t i, double val) {
    const double delta = val - vals_[i];
    vals_[i] = val;
    for (size_t k = i + 1; k < tree_.size(); k += lowbit(k)) tree_[k] += delta;
    // incremental updates accumulate rounding error; rebuilding exactly
    // every n updates keeps that bounded at amortised O(1) cost
    if (++num_updates_ > std::max<size_t>(vals_.size(), 1024)) rebuild();
  }

  // add a weight at the end
  void append(double val) {
    vals_.push_back(val);
    const size_t n = vals_.size();          // 1-based index of the new slot
    double block = val;
    for (size_t k = 1; k < lowbit(n); k <<= 1) block += tree_[n - k];
    tree_.push_back(block);
  }

  // draw an index with probability proportional to its weight
  template <typename Reng>
  size_t operator()(Reng& reng) const {
    std::uniform_real_distribution<double> unif_dist(0.0, 1.0);
    for (int attempt = 0; attempt < 100; ++attempt) {
      const size_t index = find(sum_val() * unif_dist(reng));
      // guard against landing on a zero weight through rounding error
      if (index < vals_.size() && vals_[index] > 0.0) return index;
    }
    // practically unreachable: fall back to a linear scan
    double u = sum_val() * unif_dist(reng);
    for (size_t i = 0; i < vals_.size(); ++i) {
      if (vals_[i] > 0.0) {
        if (u < vals_[i]) return i;
        u -= vals_[i];
      }
    }
    for (size_t i = vals_.size(); i-- > 0; ) if (vals_[i] > 0.0) return i;
    return 0;
  }

  // smallest 0-based index i with w_0 + ... + w_i > u
  size_t find(double u) const {
    size_t pos = 0;
    for (size_t step = highest_power_of_two(vals_.size()); step > 0; step >>= 1) {
      const size_t next = pos + step;
      if (next < tree_.size() && tree_[next] <= u) {
        pos = next;
        u -= tree_[next];
      }
    }
    return pos;  // 1-based pos + 1 -> 0-based pos
  }

private:
  std::vector<double> vals_;   // the weights
  std::vector<double> tree_;   // Fenwick tree, 1-based, tree_[0] unused
  size_t num_updates_ = 0;


  static size_t lowbit(size_t i) { return i & (~i + 1); }

  static size_t highest_power_of_two(size_t n) {
    size_t p = 1;
    while (p <= n / 2) p <<= 1;
    return n == 0 ? 0 : p;
  }

  double prefix_sum(size_t n) const {  // w_0 + ... + w_{n-1}
    double s = 0.0;
    for (size_t k = n; k > 0; k -= lowbit(k)) s += tree_[k];
    return s;
  }

  void rebuild() {  // exact O(n) construction from vals_
    tree_.assign(vals_.size() + 1, 0.0);
    for (size_t i = 1; i < tree_.size(); ++i) {
      tree_[i] += vals_[i - 1];
      const size_t parent = i + lowbit(i);
      if (parent < tree_.size()) tree_[parent] += tree_[i];
    }
    num_updates_ = 0;
  }
};
