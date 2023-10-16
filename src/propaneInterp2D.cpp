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
  bool zfix;

  int dimx = image.nrow();
  int dimy = image.ncol();

  if(zero){
    image.fill(0);
  }

  double shift = 0;

  if(FITS){
    shift = 0.5;
  }

  if(z.length() == 1){
    zfix = true;
  }else{
    zfix = false;
  }

  double sign = (type == 1) ? 1 : -1;

  for(loc = 0; loc < x.length(); loc++){
    // Rcpp::Rcout << loc << "\n";
    for (i = -1; i < 2; i++) {
      for (j = -1; j < 2; j++) {
        xloc = floor(x[loc] + i - shift);
        if(xloc >= 0 && xloc < dimx){ //since xloc will be 1 indexed
          yloc = floor(y[loc] + j - shift);
          if(yloc >= 0 && yloc < dimy){ //since yloc will be 1 indexed
            xfrac = 1 - abs(xloc - x[loc] - (shift - 0.5));
            if(xfrac >= 0 && xfrac <= 1){
              yfrac = 1 - abs(yloc - y[loc] - (shift - 0.5));
              // Rcpp::Rcout << xfrac << " " << yfrac << "\n";
              if(yfrac >= 0 && yfrac <= 1){
                if(zfix){
                  image(xloc, yloc) += z[0]*xfrac*yfrac*sign;
                }else{
                  image(xloc, yloc) += z[loc]*xfrac*yfrac*sign;
                }
              }
            }
          }
        }
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
  bool zfix;

  int dimx = image.nrow();
  int dimy = image.ncol();

  if(zero){
    image.fill(0);
  }

  double shift = 0;

  if(FITS){
    shift = 0.5;
  }

  if(z.length() == 1){
    zfix = true;
  }else{
    zfix = false;
  }

  double sign = (type == 1) ? 1 : -1;

  for(loc = 0; loc < x.length(); loc++){
    xloc = floor(x[loc] - shift);
    if(xloc >= 0 && xloc < dimx){ //since xloc will be 1 indexed
      yloc = floor(y[loc] - shift);
      if(yloc >= 0 && yloc < dimy){ //since yloc will be 1 indexed
        if(zfix){
          image(xloc, yloc) += z[0]*sign;
        }else{
          image(xloc, yloc) += z[loc]*sign;
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

  double shift = 0;

  if(FITS){
    shift = 0.5;
  }

  int sign = (type == 1) ? 1 : -1;

  for(loc = 0; loc < x.length(); loc++){
    xloc = floor(x[loc] - shift);
    if(xloc >= 0 && xloc < dimx){ //since xloc will be 1 indexed
      yloc = floor(y[loc] - shift);
      if(yloc >= 0 && yloc < dimy){ //since yloc will be 1 indexed
          image(xloc, yloc) += sign;
      }
    }
  }

  return image;
}
