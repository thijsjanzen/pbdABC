#pragma once

#include <Rcpp.h>
#include <vector>
#include <array>


inline void vector_to_numericmatrix(const std::vector< std::array< float, 4 >>& v,
                                    Rcpp::NumericMatrix& m) {
  int n_rows = v.size();
  m = Rcpp::NumericMatrix(n_rows, 4);
  for (int i = 0; i < n_rows; ++i) {
    for (int j = 0; j < 4; ++j) {
      m(i, j) = v[i][j];
    }
  }
  return;
}

inline void vector_to_numericmatrix(const std::vector< std::array< double, 4 >>& v,
                                    Rcpp::NumericMatrix& m) {
  int n_rows = v.size();
  m = Rcpp::NumericMatrix(n_rows, 4);
  for (int i = 0; i < n_rows; ++i) {
    for (int j = 0; j < 4; ++j) {
      m(i, j) = v[i][j];
    }
  }
  return;
}

inline void particle_to_numericmatrix(const std::vector< std::array<double, 10>>& v,
                                      Rcpp::NumericMatrix& m) {
  int n_rows = v.size();
  m = Rcpp::NumericMatrix(n_rows, 10);
  for (int i = 0; i < n_rows; ++i) {
    for (int j = 0; j < 10; ++j) {
      m(i, j) = v[i][j];
    }
  }
  return;
}


template <typename SPECK>
inline void update_output(std::vector< std::array<double, 10>>& out,
                   const std::vector< SPECK >& gen,
                   int iteration) {

  for (const auto& i : gen) {

    std::array<double, 10> to_add;
    to_add[0] = static_cast<double>(iteration);
    for (size_t j = 0; j < i.params_.size(); ++j) {
      to_add[j + 1] = i.params_[j];
    }
    to_add[6] = i.gamma;
    to_add[7] = i.colless;
    to_add[8] = static_cast<double>(i.num_lin);
    to_add[9] = i.weight;

    out.push_back(to_add);
  }
  return;
}

