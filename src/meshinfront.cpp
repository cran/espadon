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
// [[Rcpp::export(name = ".meshinfront")]] 
std::vector <bool> meshinfront (
    std::vector <double> pt_x,
    std::vector <double> pt_y,
    std::vector <double> pt_z,
    std::vector <double> p2_x,
    std::vector <double> p2_y,
    std::vector <double> p2_z,
    std::vector <double> u_x,
    std::vector <double> u_y,
    std::vector <double> u_z,
    std::vector <int> n_A,
    std::vector <int> n_B,
    std::vector <int> n_C
){
  // std::vector  <float> max_slice_value,
  
  int p2_idx , tr_idx, le;
  // double ABx,ABy,ABz;
  // double ACx,ACy,ACz;
  double m[9],inv_m[9], det, AP2[3], k1,k2,kd;//,k_;
  
  le = (int) p2_x.size ();
  std::vector <bool> infront(le);
  for(p2_idx = 0; p2_idx < le; p2_idx++){
    infront[p2_idx] = true;
    R_CheckUserInterrupt();
    for(tr_idx = 0; tr_idx< (int)n_A.size (); tr_idx++) {
      m[0] = pt_x[n_B[tr_idx]] -pt_x[n_A[tr_idx]];//m(0, 0)
      m[1] = pt_y[n_B[tr_idx]] -pt_y[n_A[tr_idx]];//m(1, 0)
      m[2] = pt_z[n_B[tr_idx]] -pt_z[n_A[tr_idx]];//m(2, 0)
      m[3] = pt_x[n_C[tr_idx]] -pt_x[n_A[tr_idx]];//m(0, 1)
      m[4] = pt_y[n_C[tr_idx]] -pt_y[n_A[tr_idx]];//m(1, 1)
      m[5] = pt_z[n_C[tr_idx]] -pt_z[n_A[tr_idx]];//m(2, 1)
      m[6] = -u_x[p2_idx];//m(0, 2)
      m[7] = -u_y[p2_idx];//m(1, 2)
      m[8] = -u_z[p2_idx];//m(2, 2)
      det = m[0] * (m[4] * m[8] - m[5] * m[7]) -
        m[3] * (m[1] * m[8] - m[7] * m[2]) +
        m[6] * (m[1] * m[5] - m[4] * m[2]);

      if (round(det*1.0e6)/1.0e6!=0.0) {
        inv_m[0] =  (m[4]*m[8]-m[5]*m[7])/det;
        inv_m[3] = -(m[3]*m[8]-m[6]*m[5])/det;
        inv_m[6] =  (m[3]*m[7]-m[6]*m[4])/det;
        inv_m[1] = -(m[1]*m[8]-m[7]*m[2])/det;
        inv_m[4] =  (m[0]*m[8]-m[6]*m[2])/det;
        inv_m[7] = -(m[0]*m[7]-m[1]*m[6])/det;
        inv_m[2] =  (m[1]*m[5]-m[2]*m[4])/det;
        inv_m[5] = -(m[0]*m[5]-m[2]*m[3])/det;
        inv_m[8] =  (m[0]*m[4]-m[1]*m[3])/det;

        AP2 [0] = p2_x[p2_idx] - pt_x[n_A[tr_idx]];
        AP2 [1] = p2_y[p2_idx] - pt_y[n_A[tr_idx]];
        AP2 [2] = p2_z[p2_idx] - pt_z[n_A[tr_idx]];
        k1 = round((AP2 [0] * inv_m[0] + AP2 [1] * inv_m[3] + AP2 [2] * inv_m[6])*1.0e9)/1.0e9;
        k2 = round((AP2 [0] * inv_m[1] + AP2 [1] * inv_m[4] + AP2 [2] * inv_m[7])*1.0e9)/1.0e9;
        kd = round((AP2 [0] * inv_m[2] + AP2 [1] * inv_m[5] + AP2 [2] * inv_m[8])*1.0e9)/1.0e9;
        
        if ((kd>=0.0) && (k1>=0.0) && (k2>=0.0) &&  (k1+k2<=1.0) && kd <1) {//(k_>=0.0) && 
          infront[p2_idx] = false;
          break;
        }
      }
      
    }
  }
  return(infront);
}
