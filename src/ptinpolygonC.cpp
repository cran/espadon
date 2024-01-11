#include <Rcpp.h>
using namespace Rcpp;

#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>
// #include <iostream>
#include <stdio.h>
#include <string.h>
//#include <cstdint>
#include <locale> 

int signC(double x) {
  if (x > 0.0) {
    return 1;
  } else if (x == 0.0) {
    return 0;
  } else {
    return -1;
  }
}

// [[Rcpp::export(name = ".ptinpolygonC")]] 
std::vector <int> ptinpolygonC (
    std::vector <double> point_x,
    std::vector <double> point_y, 
    std::vector <double> pol_x, 
    std::vector <double> pol_y,
    double eps = 1.0e-9){
  //il faut que point_x[0] soit le min de point_x
  unsigned int pt_idx, pol_idx, pt_le, pol_le;
  int after, s_slope, S, last_slope0, last_slope;
  double diff_x, diff1_x, diff_y,diff1_y, slope, slope_;//, sa_old;
  pt_le = (unsigned int) point_x.size ();
  pol_le = (unsigned int) pol_x.size ();
  std::vector <int> out(pt_le);
  
  last_slope0 =  signC(pol_x[pol_le - 1] - pol_x[pol_le - 2]);
  for (pt_idx = 0; pt_idx < pt_le; pt_idx++){
    R_CheckUserInterrupt();
    after = 0;
    out[pt_idx] = 0;
    // sa_old =0;
    last_slope = last_slope0;
    for (pol_idx = 0; pol_idx < (pol_le - 1); pol_idx++){
      diff_x = point_x[pt_idx] - pol_x[pol_idx];
      diff_y = point_y[pt_idx] - pol_y[pol_idx];
      diff1_y = point_y[pt_idx] - pol_y[pol_idx + 1];
      diff1_x = point_x[pt_idx] - pol_x[pol_idx + 1];
      slope_ = (pol_x[pol_idx + 1] - pol_x[pol_idx]);
      slope = slope_;
      if (fabs(slope)<=eps) slope = 0.0;
      s_slope =  signC(slope);
      if (fabs(diff_x) <= eps){

        if (fabs(diff_y) <= eps) {
          out[pt_idx] = 3; 
          break;
        } else if ((diff_y*diff1_y < 0.0) && (s_slope == 0)) { // compris entre ]y1,y2[ et pente verticale
          out[pt_idx] = 2; 
          break;
        } else if ((diff_y < 0.0) && (s_slope!=0) && (s_slope==last_slope)){
          after++;
        } 
      } else if (diff1_x * diff_x < 0.0) { // compris entre ]x1,x2[
        slope_ = (pol_y[pol_idx + 1] - pol_y[pol_idx])* diff_x - (slope_ * diff_y);
        if (fabs(slope_)<=eps) slope_ = 0.0;
        S =  signC(slope_);
        if ((S == 0) && (diff_y * diff1_y <= 0.0)){  //appartient à la droite &  compris entre [x1,x2]
          out[pt_idx] = 2; 
          break;
        } else if (s_slope == S){
          after++;
        }
      } 
      if (s_slope!=0) last_slope = s_slope;
      
    }
    if (out[pt_idx]<2) out[pt_idx]= after%2;
  }
  return(out);
}


