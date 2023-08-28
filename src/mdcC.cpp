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
  // [[Rcpp::export(name = ".mdcC")]] 
std::vector <double> mdcC (
    std::vector <double> vol3Dref, 
    std::vector <int> n_ijk, 
    std::vector <double> dxyz, 
  std::vector <int> Testi,
  std::vector <int> Testj,
  std::vector <int> Testk){
  
  int test_idx, volindex, test_le, idx;
  double dist, dist_max, di;
  test_le = (int) Testi.size ();

  std::vector <double> out(test_le);

  dist_max = sqrt(pow(n_ijk[0]*dxyz[0],2) +  pow(n_ijk[1]*dxyz[1],2)  +  pow(n_ijk[2]*dxyz[2],2));
  for (test_idx= 0; test_idx < test_le; test_idx ++){
    R_CheckUserInterrupt();
    out[test_idx] = R_NaN;
    dist = dist_max;
    for (idx=Testi[test_idx]+1;idx<n_ijk[0];idx ++){
      volindex =  idx + (Testj[test_idx] * n_ijk[0]) +  (Testk[test_idx] * n_ijk[0] * n_ijk[1]);
      di = fabs(((double)(idx-Testi[test_idx])) * dxyz[0]);
      if (di>dist) break;
      if (vol3Dref[volindex]>0) {dist=di;break;}
    }
    
    for (idx=Testi[test_idx]-1;idx>-1;idx --){
      volindex =  idx + (Testj[test_idx] * n_ijk[0]) +  (Testk[test_idx] * n_ijk[0] * n_ijk[1]);
      di = fabs(((double)(idx-Testi[test_idx])) * dxyz[0]);
      if (di>dist) break;
      if (vol3Dref[volindex]>0) {dist=di;break;}
    }
    
    for (idx=Testj[test_idx]+1;idx<n_ijk[1];idx ++){
      volindex =  Testi[test_idx] + (idx * n_ijk[0]) +  (Testk[test_idx] * n_ijk[0] * n_ijk[1]);
      di = fabs(((double)(idx-Testj[test_idx])) * dxyz[1]);
      if (di>dist) break;
      if (vol3Dref[volindex]>0) {dist=di;break;}
    }
    
    for (idx=Testj[test_idx]-1;idx>-1;idx --){
      volindex =  Testi[test_idx] + (idx * n_ijk[0]) +  (Testk[test_idx] * n_ijk[0] * n_ijk[1]);
      di = fabs(((double)(idx-Testj[test_idx])) * dxyz[1]);
      if (di>dist) break;
      if (vol3Dref[volindex]>0) {dist=di;break;}
    }
    
    for (idx=Testk[test_idx]+1;idx<n_ijk[2];idx ++){
      volindex =  Testi[test_idx] + (Testj[test_idx] * n_ijk[0]) +  (idx * n_ijk[0] * n_ijk[1]);
      di = fabs(((double)(idx-Testk[test_idx])) * dxyz[2]);
      if (di>dist) break;
      if (vol3Dref[volindex]>0) {dist=di;break;}
    }
    
    for (idx=Testk[test_idx]-1;idx>-1;idx --){
      volindex =  Testi[test_idx] + (Testj[test_idx] * n_ijk[0]) +  (idx * n_ijk[0] * n_ijk[1]);
      di = fabs(((double)(idx-Testk[test_idx])) * dxyz[2]);
      if (di>dist) break;
      if (vol3Dref[volindex]>0) {dist=di;break;}
    }
    
    if (dist!=dist_max) out[test_idx] = dist;
  }
  return(out);
}
