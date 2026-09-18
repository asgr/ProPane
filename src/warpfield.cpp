#include <Rcpp.h>
#include <vector>
using namespace Rcpp;

// Bilinear reconstruction of a warp field from a coarse lattice.
//
// The field is decomposed as   field(i, j) = affine(i, j) + residual(i, j),
// where the affine part is a least-squares fit to the coarse samples. Only the
// residual is interpolated here: subtracting the linear gradient leaves a field
// whose within-cell variation is set purely by curvature, which is what makes
// bilinear accurate to well below a pixel even on a coarse lattice.
//
// Rx/Ry are the residual samples, nc x nr in column-major order with the x
// index running fastest, i.e. laid out as matrix(v, nc, nr) over a grid built
// by rep(cxs, times = nr) / rep(cys, each = nc).
//
// cxs/cys hold 1-based lattice coordinates in R's indexing convention, matching
// the x/y values the WCS transform was evaluated at.
//
// Output is a single nx x ny x 1 x 2 array, which as.cimg() turns into the same
// 2-plane cimg that imappend() of two as.cimg()s would, so the caller skips a
// whole extra copy.
//
// Cell index and local weight are resolved by one walk per axis; the inner loop
// runs on raw pointers, with the affine term hoisted into a per-row accumulator.
// [[Rcpp::export(".warpfield_interp_cpp")]]
NumericVector warpfield_interp(NumericVector Rx, NumericVector Ry,
                               IntegerVector cxs, IntegerVector cys,
                               double a0x, double ax_x, double ax_y,
                               double a0y, double ay_x, double ay_y,
                               int nx, int ny)
{
  const int nc = cxs.size();
  const int nr = cys.size();

  if (nc < 2 || nr < 2) Rcpp::stop("coarse lattice must be at least 2x2");
  if (Rx.size() != nc * nr || Ry.size() != nc * nr)
    Rcpp::stop("residual grid does not match lattice dimensions");

  std::vector<int> ci(nx), cj(ny);
  std::vector<double> wx(nx), wy(ny);

  // Clamp the walk at the second-to-last cell so the final cell absorbs any
  // trailing pixels beyond the last interior lattice point.
  const int* cxs_p = INTEGER(cxs);
  const int* cys_p = INTEGER(cys);

  int cell = 0;
  for (int i = 0; i < nx; i++) {
    while (cell < nc - 2 && (double)(i + 1) >= (double)cxs_p[cell + 1]) cell++;
    const double x0 = (double)cxs_p[cell];
    const double x1 = (double)cxs_p[cell + 1];
    ci[i] = cell;
    double w = (x1 > x0) ? ((double)(i + 1) - x0) / (x1 - x0) : 0.0;
    wx[i] = w < 0.0 ? 0.0 : (w > 1.0 ? 1.0 : w);
  }

  cell = 0;
  for (int j = 0; j < ny; j++) {
    while (cell < nr - 2 && (double)(j + 1) >= (double)cys_p[cell + 1]) cell++;
    const double y0 = (double)cys_p[cell];
    const double y1 = (double)cys_p[cell + 1];
    cj[j] = cell;
    double w = (y1 > y0) ? ((double)(j + 1) - y0) / (y1 - y0) : 0.0;
    wy[j] = w < 0.0 ? 0.0 : (w > 1.0 ? 1.0 : w);
  }

  const double* rx = REAL(Rx);
  const double* ry = REAL(Ry);

  // One buffer holding both planes back to back: [nx*ny for x, nx*ny for y].
  // That is exactly the layout as.cimg() expects for a 2-plane cimg, so the
  // caller avoids the extra copy imappend() would make.
  const int len = nx * ny;
  NumericVector res2(len * 2);
  double* ox = REAL(res2);
  double* oy = REAL(res2) + len;

  for (int j = 0; j < ny; j++) {
    const int jb = cj[j];
    const double fy = wy[j];
    const double ly = 1.0 - fy;
    const double rowy = a0y + ay_y * (double)(j + 1);
    const double rowx = a0x + ax_y * (double)(j + 1);
    const int colbase = jb * nc;
    const int rowoff = j * nx;

    for (int i = 0; i < nx; i++) {
      const int ia = ci[i];
      const double fx = wx[i];
      const double lx = 1.0 - fx;

      const double w00 = lx * ly;
      const double w10 = fx * ly;
      const double w01 = lx * fy;
      const double w11 = fx * fy;

      const int p = ia + colbase;
      const double x = (double)(i + 1);

      ox[rowoff + i] = rowx + ax_x * x +
        (w00 * rx[p]      + w10 * rx[p + 1] +
         w01 * rx[p + nc] + w11 * rx[p + nc + 1]);
      oy[rowoff + i] = rowy + ay_x * x +
        (w00 * ry[p]      + w10 * ry[p + 1] +
         w01 * ry[p + nc] + w11 * ry[p + nc + 1]);
    }
  }

  res2.attr("dim") = IntegerVector::create(nx, ny, 1, 2);
  return res2;
}
