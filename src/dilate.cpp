#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export(".dilate_cpp")]]
IntegerMatrix dilate_cpp(IntegerMatrix segim, IntegerMatrix kern, IntegerVector expand=0){
  
  int srow = segim.nrow();
  int scol = segim.ncol();
  int krow = kern.nrow();
  int kcol = kern.ncol();
  int krow_off = ((krow - 1) / 2);
  int kcol_off = ((kcol - 1) / 2);
  int max_segim = max(segim);
  IntegerMatrix segim_new(srow, scol);
  LogicalVector seglogic  (max_segim);
  
  bool do_expand = (expand(0) > 0);

  if(do_expand){
    for(int k = 0; k < expand.length(); k++){
      if(expand(k) <= max_segim){
        seglogic(expand(k)-1) = true;
      }
    }
  }
  
  // Column-major iteration for cache locality (Rcpp matrices are column-major)
  for (int j = 0; j < scol; j++) {
    for (int i = 0; i < srow; i++) {
      if(segim(i,j) > 0){
        bool checkseg = true;
        if(do_expand){
          checkseg = seglogic(segim(i,j)-1);
        }
        if(checkseg){
          int m_start = std::max(0, krow_off - i);
          int m_end = std::min(krow, krow_off - (i - srow));
          int n_start = std::max(0, kcol_off - j);
          int n_end = std::min(kcol, kcol_off - (j - scol));
          int i_base = i - krow_off;
          int j_base = j - kcol_off;
          for (int m = m_start; m < m_end; m++) {
            for (int n = n_start; n < n_end; n++) {
              if(kern(m,n) > 0){
                if(m != krow_off || n != kcol_off){
                  int xloc = i_base + m;
                  int yloc = j_base + n;
                  if(segim(xloc,yloc) == 0){
                    if(segim_new(xloc,yloc) == 0 || segim(i,j) < segim_new(xloc,yloc)) {
                      segim_new(xloc,yloc) = segim(i,j);
                    }
                  }
                }else{
                  segim_new(i,j) = segim(i,j);
                }
              }
            }
          }
        }else{
          segim_new(i,j) = segim(i,j);
        }
      }
    }
  }
  return segim_new;
}

// IntegerMatrix dilate_cpp_old(IntegerMatrix segim, IntegerMatrix kern){
//   
//   int srow = segim.nrow();
//   int scol = segim.ncol();
//   int krow = kern.nrow();
//   int kcol = kern.ncol();
//   int krow_off = ((krow - 1) / 2);
//   int kcol_off = ((kcol - 1) / 2);
//   int maxint = std::numeric_limits<int>::max();
//   IntegerMatrix segim_new(srow, scol);
//   
//   for (int j = 0; j < scol; j++) {
//     for (int i = 0; i < srow; i++) {
//       int segID = maxint;
//       if(segim(i,j) == 0){
//         for (int n = std::max(0,kcol_off - j); n < std::min(kcol, kcol_off - (j - scol)); n++) {
//           for (int m = std::max(0,krow_off - i); m < std::min(krow, krow_off - (i - srow)); m++) {
//             if(kern(m,n) > 0){
//               int segim_segID = segim(i + m - krow_off,j + n - kcol_off);
//               if(segim_segID > 0 & segim_segID < segID) {
//                 segID = segim_segID;
//               }
//             }
//           }
//         }
//         if(segID < maxint){
//           segim_new(i,j) = segID;
//         }
//       }else{
//         segim_new(i,j) = segim(i,j);
//       }
//     }
//   }
//   
//   return segim_new;
// }
