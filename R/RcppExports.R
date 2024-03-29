# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

.dilate_cpp <- function(segim, kern, expand = 0L) {
    .Call(`_ProPane_dilate_cpp`, segim, kern, expand)
}

.mat_diff_sum <- function(mat1, mat2, scale = 1, delta_x = 0L, delta_y = 0L) {
    .Call(`_ProPane_mat_diff_sum`, mat1, mat2, scale, delta_x, delta_y)
}

.propaneInterp2D <- function(x, y, z, image, FITS = TRUE, type = 1L, zero = FALSE) {
    .Call(`_ProPane_propaneInterp2D`, x, y, z, image, FITS, type, zero)
}

.propaneBin2D <- function(x, y, z, image, FITS = TRUE, type = 1L, zero = FALSE) {
    .Call(`_ProPane_propaneBin2D`, x, y, z, image, FITS, type, zero)
}

.propaneBin2Dint <- function(x, y, image, FITS = TRUE, type = 1L, zero = FALSE) {
    .Call(`_ProPane_propaneBin2Dint`, x, y, image, FITS, type, zero)
}

.stack_image_inVar_cpp <- function(post_image, post_inVar, post_weight, pre_image, pre_inVar, pre_weight_sexp, offset, post_mask = NULL) {
    .Call(`_ProPane_stack_image_inVar`, post_image, post_inVar, post_weight, pre_image, pre_inVar, pre_weight_sexp, offset, post_mask)
}

.stack_image_cpp <- function(post_image, post_weight, pre_image, pre_weight_sexp, offset, post_mask = NULL) {
    .Call(`_ProPane_stack_image`, post_image, post_weight, pre_image, pre_weight_sexp, offset, post_mask)
}

.stack_exp_cpp <- function(post_exp, pre_exp, offset) {
    .Call(`_ProPane_stack_exp`, post_exp, pre_exp, offset)
}

.stack_exp_mask_cpp <- function(post_exp, pre_exp, offset, post_mask) {
    .Call(`_ProPane_stack_exp_mask`, post_exp, pre_exp, offset, post_mask)
}

