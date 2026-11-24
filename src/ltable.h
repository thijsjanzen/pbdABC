#pragma once
#include <algorithm>
#include <array>
#include <cmath>
#include <numeric>
#include <unordered_map>
#include <vector>


using ltable = std::vector< std::array<double, 4>>;

enum species_property {birth_time, parent, id, death_time};



inline void sort_by_time(ltable& ltab) {
  std::sort(ltab.begin(), ltab.end(), [&](const auto& a, const auto& b) {
    return a[species_property::birth_time] > b[species_property::birth_time];
  });
}

inline std::vector<double> brts_from_ltable(ltable L) {
  sort_by_time(L);

  std::vector<double> out;
  out.reserve(L.size() - 1);
  for (size_t i = 1; i < L.size(); ++i) {
      out.push_back(L[i][species_property::birth_time]);
  }

  return out;
}


// Reconstructed tree: drop extinct lineages from an L table.
//
// When a lineage goes extinct but has surviving descendants, the descendant
// branch that split off from it last continues its identity (after that last
// split the lineage's own branch dies, so the node disappears from the
// reconstructed tree). So every extant lineage carries the identity (id,
// birth time and parent) of the oldest ancestor it continues.
//
// Single pass from youngest to oldest lineage plus one sort: O(n log n).
// (The previous version re-sorted the table once per replaced row, which was
// O(n^2 log n).) Rows have columns birth_time (age), parent, id, death_time
// (-1 if extant). Returns an empty table if fewer than two lineages survive.
inline ltable drop_extinct(const ltable& L) {
  using sp = species_property;
  const size_t n = L.size();
  const size_t none = n;

  // row of each id
  std::unordered_map<long, size_t> row_of;
  row_of.reserve(2 * n);
  for (size_t i = 0; i < n; ++i) row_of[std::lround(L[i][sp::id])] = i;
  auto parent_row = [&](size_t i) {
    auto it = row_of.find(std::lround(L[i][sp::parent]));
    return it == row_of.end() ? none : it->second;
  };

  // youngest first (smallest age), so all children of a lineage are visited
  // before the lineage itself; ties (the two crown lineages) are broken by
  // row, later rows first, since children are always appended after parents
  std::vector<size_t> order(n);
  std::iota(order.begin(), order.end(), 0);
  std::sort(order.begin(), order.end(), [&](size_t a, size_t b) {
    if (L[a][sp::birth_time] != L[b][sp::birth_time]) {
      return L[a][sp::birth_time] < L[b][sp::birth_time];
    }
    return a > b;
  });

  // carrier[i]: the extant row that continues lineage i (none if lineage i
  // has no surviving descendants); latest_child[i]: its youngest child with
  // surviving descendants
  std::vector<size_t> carrier(n, none);
  std::vector<size_t> latest_child(n, none);
  for (size_t i : order) {
    if (L[i][sp::death_time] < 0) {
      carrier[i] = i;
    } else if (latest_child[i] != none) {
      carrier[i] = carrier[latest_child[i]];
    }
    if (carrier[i] == none) continue;
    const size_t p = parent_row(i);
    if (p != none) {
      size_t& lc = latest_child[p];
      if (lc == none || L[i][sp::birth_time] < L[lc][sp::birth_time]) lc = i;
    }
  }

  // for each extant row: the oldest ancestor whose identity it carries
  std::vector<size_t> top(n, none);
  for (size_t y = 0; y < n; ++y) {
    if (carrier[y] != y) continue;
    size_t a = y;
    for (size_t p = parent_row(a); p != none && carrier[p] == y; p = parent_row(a)) {
      a = p;
    }
    top[y] = a;
  }

  // one row per extant lineage, with the identity it carries
  struct entry { std::array<double, 4> row; size_t top_row; };
  std::vector<entry> rows;
  for (size_t y = 0; y < n; ++y) {
    if (top[y] == none) continue;
    const size_t a = top[y];
    const size_t p = parent_row(a);
    // the parent lineage is continued by carrier[p], under that row's identity
    const double parent_id = p == none ? 0.0 : L[top[carrier[p]]][sp::id];
    rows.push_back({{L[a][sp::birth_time], parent_id, L[a][sp::id],
                     L[y][sp::death_time]}, a});
  }
  if (rows.size() < 2) return ltable();

  // oldest first; the two crown lineages (equal age) in their original order
  std::sort(rows.begin(), rows.end(), [](const entry& a, const entry& b) {
    if (a.row[sp::birth_time] != b.row[sp::birth_time]) {
      return a.row[sp::birth_time] > b.row[sp::birth_time];
    }
    return a.top_row < b.top_row;
  });

  // renumber ids 1..m in order of birth, keeping the sign (crown side)
  std::unordered_map<long, double> new_id;
  new_id.reserve(2 * rows.size());
  ltable out;
  out.reserve(rows.size());
  for (size_t i = 0; i < rows.size(); ++i) {
    const double old_id = rows[i].row[sp::id];
    const double id_i = old_id < 0 ? -static_cast<double>(i + 1)
                                   : static_cast<double>(i + 1);
    new_id[std::lround(old_id)] = id_i;
    out.push_back(rows[i].row);
    out.back()[sp::id] = id_i;
  }
  for (auto& r : out) {
    auto it = new_id.find(std::lround(r[sp::parent]));
    r[sp::parent] = it == new_id.end() ? 0.0 : it->second;
  }

  out[1][sp::parent] = out[0][sp::id];
  out[0][sp::parent] = 0;

  return out;
}
