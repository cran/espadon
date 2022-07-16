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
// [[Rcpp::export(name = ".getvaluefromijkC")]] 
std::vector <double> getvaluefromijkC (
    std::vector <double> vol3D,
    bool interpolate,
    std::vector <double> i,
    std::vector <double> j,
    std::vector <double> k,
    std::vector  <int> k_idx,
    std::vector  <int> k_loc,
    std::vector  <int> n_ijk){
  
  long idx, le, le_k, le_map;
  int round_i, round_j, round_k;
  bool test;
  // int test_idx;
  double pi,pj,pk;
  
  le = (long) i.size ();
  le_k = (long) k_idx.size ();
  le_map = n_ijk[0] * n_ijk[1];
  
  std::vector <double> value(le);
  for (idx= 0; idx < le; idx ++) value[idx] = R_NaN;
  
  for (idx= 0; idx < le; idx ++){
    
    
    
    if (!interpolate){
      if (i[idx]<0) {
        round_i = (int) (i[idx] -  0.5);
      } else { 
        round_i = (int) (i[idx] +  0.5);
      }
      if (j[idx]<0) {
        round_j = (int) (j[idx] -  0.5);
      } else  {
        round_j = (int) (j[idx] +  0.5);
      }
      if (k[idx]<0) {
        round_k = (int) (k[idx] -  0.5);
      } else  {
        round_k = (int) (k[idx] +  0.5);
      }
      // Rcout<<std::to_string(round_i)<<" "<<std::to_string(round_j)<<" "<< std::to_string(round_k) <<  "\n";
      test = (round_k>=0) && (round_k<=k_idx[le_k-1]);
      
      if (test) 
        test = (k_idx[round_k]==round_k) && (round_i>=0) &&  (round_i<n_ijk[0]) && (round_j>=0) && (round_j<n_ijk[1]);
      
      if (test) 
        value[idx] = vol3D[(round_i + (round_j * n_ijk[0]) +  (k_loc[round_k] * le_map))];
      
      
    } else {
      if (i[idx]<0) {
        round_i = (int) i[idx] - 1;
      } else { 
        round_i = (int) i[idx];
      }
      if (j[idx]<0) {
        round_j = (int) j[idx] - 1;
      } else  {
        round_j = (int) j[idx];
      }
      if (k[idx]<0) {
        round_k = (int) k[idx] - 1;
      } else  {
        round_k = (int) k[idx];
      }
      
      // Rcout<<std::to_string(round_i)<<" "<<std::to_string(round_j)<<" "<<std::to_string(round_k)<<" "<< k[idx]<<" ";
      test = (round_k>=0) && (round_k<k_idx[le_k-1]);
      if (test) test = (k_idx[round_k]==round_k) && (k_idx[round_k+1]==round_k+1) && 
        (round_i>=0) && (round_i+1<n_ijk[0]) && (round_j>=0) && (round_j+1<n_ijk[1]);
      
      if (test) {
        // RCout<<std::to_string((int)(-1.5))<<" "<<std::to_string((int)1.9)<<" ";
        
        pi = i[idx]-round_i;
        pj = j[idx]-round_j;
        pk = k[idx]-round_k;
        
        round_k = k_loc[round_k];
        
        // Rcout<<std::to_string(pi)<<" "<<std::to_string(pj)<<" "<<std::to_string(pk);
        value[idx] = 
          vol3D[(round_i + (round_j * n_ijk[0]) +  (round_k * le_map))] * (1-pi) * (1-pj) * (1-pk) +
          vol3D[(round_i + 1 + (round_j * n_ijk[0]) +  (round_k * le_map))]* pi * (1-pj) * (1-pk) +
          vol3D[(round_i + ((round_j+1) * n_ijk[0]) +  (round_k * le_map))] * (1-pi) * pj * (1-pk) +
          vol3D[(round_i + 1 + ((round_j+1) * n_ijk[0]) +  (round_k * le_map))] * pi * pj * (1-pk) +
          
          vol3D[(round_i + (round_j * n_ijk[0]) +  ((round_k+1) * le_map))] * (1-pi) * (1-pj) * pk +
          vol3D[(round_i + 1 + (round_j * n_ijk[0]) +  ((round_k+1) * le_map))]* pi * (1-pj) * pk +
          vol3D[(round_i + ((round_j+1) * n_ijk[0]) +  ((round_k+1) * le_map))] * (1-pi) * pj * pk +
          vol3D[(round_i + 1 + ((round_j+1) * n_ijk[0]) +  ((round_k+1) * le_map))] * pi * pj * pk;
        
      }
      // Rcout<<"\n";
      
    }
    
  }
  
  
  return (value);
}
