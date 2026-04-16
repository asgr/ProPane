#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export(".propaneInterp2D")]]
NumericMatrix propaneInterp2D(NumericVector x, NumericVector y, NumericVector z,
                              NumericMatrix image,
                              bool FITS=true, int type=1L, bool zero=false) {

  int i, j, loc, xloc, yloc;
  double xfrac, yfrac;

  int dimx = image.nrow();
  int dimy = image.ncol();

  if(zero){
    image.fill(0);
  }

  double shift = FITS ? 0.5 : 0.0;
  double shift_adj = shift - 0.5;
  bool zfix = (z.length() == 1);
  double sign = (type == 1) ? 1.0 : -1.0;

  for(loc = 0; loc < x.length(); loc++){
    double zval = zfix ? z[0] : z[loc];
    double zval_sign = zval * sign;
    // Precompute the base floor position once per point (not 9 times)
    double xbase = x[loc] - shift;
    double ybase = y[loc] - shift;
    int xfloor_base = (int)floor(xbase);
    int yfloor_base = (int)floor(ybase);

    for (i = -1; i < 2; i++) {
      xloc = xfloor_base + i;
      if(xloc < 0 || xloc >= dimx) continue;
      xfrac = 1.0 - fabs(xloc - x[loc] - shift_adj);
      if(xfrac <= 0.0 || xfrac > 1.0) continue;
      for (j = -1; j < 2; j++) {
        yloc = yfloor_base + j;
        if(yloc < 0 || yloc >= dimy) continue;
        yfrac = 1.0 - fabs(yloc - y[loc] - shift_adj);
        if(yfrac <= 0.0 || yfrac > 1.0) continue;
        image(xloc, yloc) += zval_sign * xfrac * yfrac;
      }
    }
  }

  return image;
}

// [[Rcpp::export(".propaneBin2D")]]
NumericMatrix propaneBin2D(NumericVector x, NumericVector y, NumericVector z,
                              NumericMatrix image,
                              bool FITS=true, int type=1L, bool zero=false) {

  int loc, xloc, yloc;

  int dimx = image.nrow();
  int dimy = image.ncol();

  if(zero){
    image.fill(0);
  }

  double shift = FITS ? 0.5 : 0.0;
  bool zfix = (z.length() == 1);
  double sign = (type == 1) ? 1.0 : -1.0;

  if(zfix){
    double zval_sign = z[0] * sign;
    for(loc = 0; loc < x.length(); loc++){
      xloc = (int)floor(x[loc] - shift);
      if(xloc >= 0 && xloc < dimx){
        yloc = (int)floor(y[loc] - shift);
        if(yloc >= 0 && yloc < dimy){
          image(xloc, yloc) += zval_sign;
        }
      }
    }
  }else{
    for(loc = 0; loc < x.length(); loc++){
      xloc = (int)floor(x[loc] - shift);
      if(xloc >= 0 && xloc < dimx){
        yloc = (int)floor(y[loc] - shift);
        if(yloc >= 0 && yloc < dimy){
          image(xloc, yloc) += z[loc] * sign;
        }
      }
    }
  }

  return image;
}

// [[Rcpp::export(".propaneBin2Dint")]]
IntegerMatrix propaneBin2Dint(NumericVector x, NumericVector y, IntegerMatrix image,
                           bool FITS=true, int type=1L, bool zero=false) {

  int loc, xloc, yloc;

  int dimx = image.nrow();
  int dimy = image.ncol();

  if(zero){
    image.fill(0);
  }

  double shift = FITS ? 0.5 : 0.0;
  int sign = (type == 1) ? 1 : -1;

  for(loc = 0; loc < x.length(); loc++){
    xloc = (int)floor(x[loc] - shift);
    if(xloc >= 0 && xloc < dimx){
      yloc = (int)floor(y[loc] - shift);
      if(yloc >= 0 && yloc < dimy){
          image(xloc, yloc) += sign;
      }
    }
  }

  return image;
}
