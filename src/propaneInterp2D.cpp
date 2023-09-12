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
  float xfrac, yfrac;
  bool zfix;

  int dimx = image.nrow();
  int dimy = image.ncol();

  if(zero){
    for (i = 0; i < dimx; i++) {
      for (j = 0; j < dimy; j++) {
        image(i, j) = 0;
      }
    }
  }

  float shift = 0;

  if(FITS){
    shift = 0.5;
  }

  if(z.length() == 1){
    zfix = true;
  }else{
    zfix = false;
  }


  for(loc = 0; loc < x.length(); loc++){
    // Rcpp::Rcout << loc << "\n";
    for (i = -1; i < 2; i++) {
      for (j = -1; j < 2; j++) {
        xloc = ceil(x[loc] + i - shift);
        if(xloc >= 1 & xloc <= dimx){ //since xloc will be 1 indexed
          yloc = ceil(y[loc] + j - shift);
          if(yloc >= 1 & yloc <= dimy){ //since yloc will be 1 indexed
            xfrac = 1 - abs(xloc - x[loc] - (0.5 - shift));
            if(xfrac >= 0 & xfrac <= 1){
              yfrac = 1 - abs(yloc - y[loc] - (0.5 - shift));
              // Rcpp::Rcout << xfrac << " " << yfrac << "\n";
              if(yfrac >= 0 & yfrac <= 1){
                if(type==1){
                  // Rcpp::Rcout << z[loc]*xfrac*yfrac << "\n";
                  if(zfix){
                    image(xloc - 1, yloc - 1) += z[0]*xfrac*yfrac;
                  }else{
                    image(xloc - 1, yloc - 1) += z[loc]*xfrac*yfrac;
                  }
                }
                if(type==2){
                  if(zfix){
                    image(xloc - 1, yloc - 1) -= z[0]*xfrac*yfrac;
                  }else{
                    image(xloc - 1, yloc - 1) -= z[loc]*xfrac*yfrac;
                  }
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

  int i, j, loc, xloc, yloc;
  float xfrac, yfrac;
  bool zfix;

  int dimx = image.nrow();
  int dimy = image.ncol();

  if(zero){
    for (i = 0; i < dimx; i++) {
      for (j = 0; j < dimy; j++) {
        image(i, j) = 0;
      }
    }
  }

  float shift = 0;

  if(FITS){
    shift = 0.5;
  }

  if(z.length() == 1){
    zfix = true;
  }else{
    zfix = false;
  }


  for(loc = 0; loc < x.length(); loc++){
    xloc = ceil(x[loc] - shift);
    if(xloc >= 1 & xloc <= dimx){ //since xloc will be 1 indexed
      yloc = ceil(y[loc] - shift);
      if(yloc >= 1 & yloc <= dimy){ //since yloc will be 1 indexed
        if(zfix){
          image(xloc - 1, yloc - 1) += z[0];
        }else{
          image(xloc - 1, yloc - 1) += z[loc];
        }
      }
    }
  }

  return image;
}


// for(ix in -1:1){
//   for(jy in -1:1){
//     xfrac = (1 - abs((xysub_loc[,1] + ix) - image_pre_fix_xysub[,1]))
//     yfrac = (1 - abs((xysub_loc[,2] + jy) - image_pre_fix_xysub[,2]))
//
//     if(ix != 0L | jy != 0L){ #since [0,0] is the centre pixel, which will always have some flux share between 0-1
//       xfrac[xfrac < 0] = 0
//       xfrac[xfrac > 1] = 1
//       yfrac[yfrac < 0] = 0
//       yfrac[yfrac > 1] = 1
//     }
//
//     subset = cbind(xysub_loc[,1] + ix, xysub_loc[,2] + jy)
//       subset_dup = duplicated(subset)
//
//       if(any(subset_dup)){
//         costmat[subset[!subset_dup,]] = costmat[subset[!subset_dup,]] - (flux[!subset_dup]*xfrac[!subset_dup]*yfrac[!subset_dup])
//         for(k in which(subset_dup)){
//           costmat[subset[k,]] = costmat[subset[k,]] - (flux[k]*xfrac[k]*yfrac[k])
//         }
//       }else{
//         costmat[subset] = costmat[subset] - (flux*xfrac*yfrac)
//       }
//   }
// }
