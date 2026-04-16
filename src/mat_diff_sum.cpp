#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export(.mat_diff_sum)]]
double mat_diff_sum(Rcpp::NumericMatrix mat1, Rcpp::NumericMatrix mat2, double scale = 1, int delta_x = 0, int delta_y = 0){
  int nrow = mat1.nrow();
  int ncol = mat1.ncol();
  double mat_sum = 0;
  double scale_inv = 1.0 / scale;

  // Precompute valid loop bounds to avoid per-element boundary checks
  int i_start = std::max(0, delta_x);
  int i_end = std::min(nrow, nrow + delta_x);
  int j_start = std::max(0, delta_y);
  int j_end = std::min(ncol, ncol + delta_y);

  // Column-major iteration for cache locality (Rcpp matrices are column-major)
  for (int j = j_start; j < j_end; j++) {
    for (int i = i_start; i < i_end; i++) {
      double v1 = mat1(i, j);
      double v2 = mat2(i - delta_x, j - delta_y);
      if(std::isfinite(v1) && std::isfinite(v2)){
        double diff = (v1 - v2) * scale_inv;
        mat_sum += diff * diff;
      }
    }
  }
  return mat_sum;
}
