#include <Rcpp.h>
using namespace Rcpp;

template<typename T, typename U>
static void check_same_size(T a, U b, const char *name)
{
  if (a.nrow() != b.nrow() || a.ncol() != b.ncol())
    Rcpp::stop(
      "%s has non-conforming size: (%d, %d) != (%d, %d))",
      name, a.nrow(), a.ncol(), b.nrow(), b.ncol()
    );
}

template<typename T>
static void check_same_size(T a, Nullable<LogicalMatrix> b, const char *name)
{
  if (b.isNull())
    return;
  check_same_size(a, LogicalMatrix(b), name);
}

static void check_required_size(NumericMatrix pre, IntegerVector offset, NumericMatrix post)
{
  if (offset.size() < 2)
    Rcpp::stop("offset has size < 2");
  auto offset_x = offset[0] - 1;
  auto offset_y = offset[1] - 1;
  if (post.nrow() < pre.nrow() + offset_x || post.ncol() < pre.ncol() + offset_y)
    Rcpp::stop(
      "post image has smaller size than required: (%d, %d) < (%d, %d) + [%d, %d])",
      post.nrow(), post.ncol(), pre.nrow(), pre.ncol(), offset_x, offset_y
    );
}

static bool is_mask_set(Nullable<LogicalMatrix> mask, int i, int j)
{
  if (mask.isNull())
    return false;
  return LogicalMatrix(mask)(i, j);
}

template<typename WeightGetter>
static SEXP _stack_image_inVar(NumericMatrix post_image, NumericMatrix post_inVar, IntegerMatrix post_weight,
                               NumericMatrix pre_image, NumericMatrix pre_inVar, WeightGetter weight, IntegerVector offset,
                               Nullable<LogicalMatrix> post_mask)
{
  check_same_size(post_image, post_inVar, "post_inVar");
  check_same_size(post_image, post_weight, "post_weight");
  check_same_size(post_image, post_mask, "post_mask");
  check_same_size(pre_image, pre_inVar, "pre_inVar");
  check_required_size(pre_image, offset, post_image);
  int post_j = offset[1] - 1;
  for (int j = 0; j < pre_image.ncol(); j++, post_j++) {
    int post_i = offset[0] - 1;
    for (int i = 0; i < pre_image.nrow(); i++, post_i++) {
      auto pre_image_pix = pre_image(i,j);
      auto pre_inVar_pix = pre_inVar(i,j);
      auto do_stack = R_finite(pre_image_pix + pre_inVar_pix) && (pre_inVar_pix > 0) && !is_mask_set(post_mask, post_i, post_j);
      if (do_stack) {
        post_image(post_i, post_j) += pre_image_pix * pre_inVar_pix;
        post_inVar(post_i, post_j) += pre_inVar_pix;
        post_weight(post_i, post_j) += weight(i,j);
      }
    }
  }
  return R_NilValue;
}



template<typename WeightGetter>
static SEXP _stack_image(NumericMatrix post_image, IntegerMatrix post_weight,
                         NumericMatrix pre_image, WeightGetter weight, IntegerVector offset,
                         Nullable<LogicalMatrix> post_mask)
{
  check_same_size(post_image, post_weight, "post_weight");
  check_same_size(post_image, post_mask, "post_mask");
  check_required_size(pre_image, offset, post_image);
  int post_j = offset[1] - 1;
  for (int j = 0; j < pre_image.ncol(); j++, post_j++) {
    int post_i = offset[0] - 1;
    for (int i = 0; i < pre_image.nrow(); i++, post_i++) {
      auto pre_image_pix = pre_image(i,j);
      auto do_stack = R_finite(pre_image_pix) && !is_mask_set(post_mask, post_i, post_j);
      if (do_stack) {
        auto weight_pix = weight(i,j);
        post_image(post_i, post_j) += pre_image_pix * weight_pix;
        post_weight(post_i, post_j) += weight_pix;
      }
    }
  }
  return R_NilValue;
}

// [[Rcpp::export(".stack_image_inVar_cpp")]]
SEXP stack_image_inVar(NumericMatrix post_image, NumericMatrix post_inVar, IntegerMatrix post_weight,
                       NumericMatrix pre_image, NumericMatrix pre_inVar, SEXP pre_weight_sexp, IntegerVector offset,
                       Nullable<LogicalMatrix> post_mask = R_NilValue)
{
  if (Rf_isMatrix(pre_weight_sexp)) {
    NumericMatrix pre_weight(pre_weight_sexp);
    check_same_size(pre_image, pre_weight, "pre_weight");
    return _stack_image_inVar(post_image, post_inVar, post_weight, pre_image, pre_inVar, [&](int i, int j) { return pre_weight(i, j); }, offset, post_mask);
  }
  int pre_weight = Rf_asInteger(pre_weight_sexp);
  return _stack_image_inVar(post_image, post_inVar, post_weight, pre_image, pre_inVar, [&](int, int) { return pre_weight; }, offset, post_mask);
}

// [[Rcpp::export(".stack_image_cpp")]]
SEXP stack_image(NumericMatrix post_image, IntegerMatrix post_weight,
                 NumericMatrix pre_image, SEXP pre_weight_sexp, IntegerVector offset,
                 Nullable<LogicalMatrix> post_mask = R_NilValue)
{
  if (Rf_isMatrix(pre_weight_sexp)) {
    NumericMatrix pre_weight(pre_weight_sexp);
    check_same_size(pre_image, pre_weight, "pre_weight");
    return _stack_image(post_image, post_weight, pre_image, [&](int i, int j) { return pre_weight(i, j); }, offset, post_mask);
  }
  int pre_weight = Rf_asInteger(pre_weight_sexp);
  return _stack_image(post_image, post_weight, pre_image, [&](int, int) { return pre_weight; }, offset, post_mask);
}

// [[Rcpp::export(".stack_exp_cpp")]]
SEXP stack_exp(NumericMatrix post_exp, NumericMatrix pre_exp, IntegerVector offset){
  int do_stack;

  for (int j = 0; j < pre_exp.ncol(); j++) {
    for (int i = 0; i < pre_exp.nrow(); i++) {
      do_stack = R_finite(pre_exp(i,j));
      if(do_stack){
        post_exp(i + offset[0] - 1, j + offset[1] - 1) += pre_exp(i,j);
      }
    }
  }
  return R_NilValue;
}

// [[Rcpp::export(".stack_exp_mask_cpp")]]
SEXP stack_exp_mask(NumericMatrix post_exp, NumericMatrix pre_exp, IntegerVector offset,
                    LogicalMatrix post_mask){
  int do_stack;

  for (int j = 0; j < pre_exp.ncol(); j++) {
    for (int i = 0; i < pre_exp.nrow(); i++) {
      do_stack = R_finite(pre_exp(i,j)) && (post_mask(i + offset[0] - 1, j + offset[1] - 1) == false);
      if(do_stack){
        post_exp(i + offset[0] - 1, j + offset[1] - 1) += pre_exp(i,j);
      }
    }
  }
  return R_NilValue;
}

