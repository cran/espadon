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
  long round_i, round_j, round_k;
  // long kn;
  bool testi,testj,testk;
  // int test_idx;
  double pi,pj,pk;
  
  le = (long) i.size ();
  le_k = (long) k_idx.size ();
  le_map = n_ijk[0] * n_ijk[1];
  
  std::vector <double> value(le);
  std::vector <double> d(7);
  // for (idx= 0; idx < le; idx ++) value[idx] = R_NaN;
  
  for (idx= 0; idx < le; idx ++){
    value[idx] = R_NaN;
    
    round_i = (int) (i[idx] +  0.5);
    round_j = (int) (j[idx] +  0.5);
    round_k = (int) (k[idx] +  0.5);
    // Rcout<<std::to_string(round_i)<<" "<<std::to_string(round_j)<<" "<<std::to_string(round_k)<<"\n";
    if (!interpolate){
      if ((round_i >= 0) && (round_j >= 0) && (round_k >= k_idx[0]) && 
          (round_i < n_ijk[0]) && (round_j < n_ijk[1]) && 
          (round_k <= k_idx[le_k-1])){
        value[idx] = vol3D[(round_i + (round_j * n_ijk[0]) +  (k_loc[round_k] * le_map))];
      }
    } else {
      round_i = (int) i[idx];
      round_j = (int) j[idx];
      round_k = (int) k[idx];
      
      pi = ((int)((i[idx]-round_i)*1.0e6))/1.0e6;
      pj = ((int)((j[idx]-round_j)*1.0e6))/1.0e6;
      pk = ((int)((k[idx]-round_k)*1.0e6))/1.0e6;
      // Rcout<<std::to_string(round_i)<<" "<<std::to_string(pi)<<" "<<std::to_string(i[idx])<<"\n";
      testi = ((pi>0) && (round_i >= 0) && (round_i< n_ijk[0]-1)) || 
        ((pi==0) && (round_i >= 0) && (round_i< n_ijk[0]));
      testj = ((pj>0) && (round_j >= 0) && (round_j< n_ijk[1]-1)) || 
        ((pj==0) && (round_j >= 0) && (round_j< n_ijk[1]));
      testk = ((pk>0) && (round_k >= k_idx[0]) && (round_k< k_idx[le_k-1])) || 
        ((pk==0) && (round_k >= k_idx[0]) && (round_k<= k_idx[le_k-1]));
      
      if (testi && testj && testk){

        // test =  ((pk==0) && (round_k<=k_idx[le_k-1])) || (round_k<k_idx[le_k-1]);
        
       
        // if (round_k<k_idx[le_k-1]) {kn = k_idx[round_k+1];} else {kn = le_k;}
        // if (test) 
        //   test = (k_idx[round_k]==round_k) && ((pk==0) ||  (kn == round_k+1)) &&
        //     (((pi==0) && (round_i<n_ijk[0])) || (round_i+1<n_ijk[0])) && 
        //     (((pj==0) && (round_j<n_ijk[1])) || (round_j+1<n_ijk[1]));
        

        // if (test) {
          
          value[idx] = vol3D[(round_i + (round_j * n_ijk[0]) +  (k_loc[round_k] * le_map))] * (1-pi) * (1-pj) * (1-pk);
          if (pi!=0) value[idx] = value[idx] + vol3D[(round_i + 1 + (round_j * n_ijk[0]) +  (k_loc[round_k] * le_map))]* pi * (1-pj) * (1-pk);
          if (pj!=0) value[idx] = value[idx] + vol3D[(round_i + ((round_j+1) * n_ijk[0]) +  (k_loc[round_k] * le_map))] * (1-pi) * pj * (1-pk);
          if ((pi!=0) && (pj!=0)) value[idx] = value[idx] + vol3D[(round_i + 1 + ((round_j+1) * n_ijk[0]) +  (k_loc[round_k] * le_map))] * pi * pj * (1-pk);
          if (pk!=0){
            value[idx] = value[idx] + vol3D[(round_i + (round_j * n_ijk[0]) +  ((k_loc[round_k]+1) * le_map))] * (1-pi) * (1-pj) * pk;
            if (pi!=0) value[idx] = value[idx] + vol3D[(round_i + 1 + (round_j * n_ijk[0]) +  ((k_loc[round_k]+1) * le_map))]* pi * (1-pj) * pk;
            if (pj!=0) value[idx] = value[idx] + vol3D[(round_i + ((round_j+1) * n_ijk[0]) +  ((k_loc[round_k]+1) * le_map))] * (1-pi) * pj * pk;
            if ((pi!=0) && (pj!=0)) value[idx] = value[idx] + vol3D[(round_i + 1 + ((round_j+1) * n_ijk[0]) +  ((k_loc[round_k]+1) * le_map))] * pi * pj * pk;
          }
          
        // }
      }
    }
  }
  return (value);
}
