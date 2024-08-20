#pragma once
#include <vector>
#include <array>


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


inline ltable drop_extinct(const ltable& L) {
  ltable out;
  out.reserve(L.size() - 2);
  for (size_t i = 0; i < L.size(); ++i) {
    if (L[i][species_property::death_time] < 0) {
      out.push_back(L[i]);
    }
  }

  // sort by branching time:
  sort_by_time(out);

  // now we have a broken Ltable, we need to fill in parents of parents
  for (int i = out.size() - 1; i >= 2; ) {
    auto parent = out[i][species_property::parent];

    if (parent == 0) break; // we are at the root.

    bool parent_found = false;

    for (size_t j = 0; j < i; ++j) {
      if (out[j][species_property::id] == parent) {
        parent_found = true;
        break;
      }
    }
    if (!parent_found) {
      size_t j = 0;
      for (; j < L.size(); ++j) {
        if (L[j][species_property::id] == parent) {
          break;
        }
      }
      auto old_id = out[i][species_property::id];

      out[i][species_property::parent]     = L[j][species_property::parent];
      out[i][species_property::birth_time] = L[j][species_property::birth_time];
      out[i][species_property::id]         = L[j][species_property::id];

      auto new_id = out[i][species_property::id];

      for (size_t k = i; k < out.size(); ++k) {
        if (out[k][species_property::parent] == old_id)
          out[k][species_property::parent] = new_id;
      }
    } else {
      --i;
    }

    sort_by_time(out);
  }

  out[0] = L[0];
  out[1] = L[1]; // to avoid sorting effects


  // and now we need to renumber to ensure consistent numbering
  for (int i = 0; i < out.size(); ++i) {
    auto id = std::abs(out[i][species_property::id]);
    auto old_id = out[i][species_property::id];
    if (id != (i + 1)) {
      int new_id = i + 1;
      if (out[i][species_property::id] < 0) new_id *= -1;

      out[i][species_property::id] = new_id;

      for (size_t j = i; j < out.size(); ++j) {
        if (out[j][species_property::parent] == old_id) {
          out[j][species_property::parent] = new_id;
        }
      }
    }
  }

  return out;
}
