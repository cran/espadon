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
    std::vector  <int> n_ijk,
    std::vector <double> s_ijk){
  long idx, le, le_k, le_map;
  bool testi,testj,testk,fill_i,fill_j,fill_k;
  std::vector <long> round_i(2),round_j(2),round_k(2);
  std::vector <double> rg_i(2),rg_j(2),rg_k(2);
  std::vector <double> pi(2),pj(2),pk(2);
  
  le = (long) i.size ();
  le_k = (long) k_idx.size ();
  le_map = n_ijk[0] * n_ijk[1];
  
  std::vector <double> value(le);
  
  for (idx= 0; idx < le; idx ++){
    value[idx] = R_NaN;
    
    round_i[0] = (int) (i[idx] +  0.5);
    round_j[0] = (int) (j[idx] +  0.5);
    round_k[0] = (int) (k[idx] +  0.5);
    // Rcout<<std::to_string(round_i)<<" "<<std::to_string(round_j)<<" "<<std::to_string(round_k)<<"\n";
    if (!interpolate){
      if ((round_i[0] >= 0) && (round_j[0] >= 0) && (round_k[0] >= k_idx[0]) && 
          (round_i[0] < n_ijk[0]) && (round_j[0] < n_ijk[1]) && 
          (round_k[0] <= k_idx[le_k-1])){
        value[idx] = vol3D[(round_i[0] + (round_j[0] * n_ijk[0]) +  (k_loc[round_k[0]] * le_map))];
      }
    } else {
      rg_i[0] = (double)i[idx] - s_ijk[0]/2.0; rg_i[1] = (double)i[idx] + s_ijk[0]/2.0;
      rg_j[0] = (double)j[idx] - s_ijk[1]/2.0; rg_j[1] = (double)j[idx] + s_ijk[1]/2.0;
      rg_k[0] = (double)k[idx] - s_ijk[2]/2.0; rg_k[1] = (double)k[idx] + s_ijk[2]/2.0;
      
      round_i[0] = (int) (rg_i[0] +  0.5); round_i[1] = (int) (rg_i[1] +  0.5);
      round_j[0] = (int) (rg_j[0] +  0.5); round_j[1] = (int) (rg_j[1] +  0.5);
      round_k[0] = (int) (rg_k[0] +  0.5); round_k[1] = (int) (rg_k[1] +  0.5);
      
      pi[0] = (double)round_i[0] + 0.5 - rg_i[0]; pi[1] = rg_i[1] - (double)round_i[1] + 0.5;
      pj[0] = (double)round_j[0] + 0.5 - rg_j[0]; pj[1] = rg_j[1] - (double)round_j[1] + 0.5;
      pk[0] = (double)round_k[0] + 0.5 - rg_k[0]; pk[1] = rg_k[1] - (double)round_k[1] + 0.5;
      
      testi = (rg_i[0] >= -0.5) && (rg_i[1]<= n_ijk[0]-0.5);
      testj = (rg_j[0] >= -0.5) && (rg_j[1]<= n_ijk[1]-0.5);
      testk = (rg_k[0] >= k_idx[0]-0.5) && (rg_k[1]<= k_idx[le_k - 1] + 0.5);
      

      // testi = (round_i[0] >= 0) && (((pi[0] < 1) && (round_i[1] < n_ijk[0])) || 
      //   ((pi[0] == 1) && (round_i[0] < n_ijk[0])));
      // testj = (round_j[0] >= 0) && (((pj[0] < 1) && (round_j[1] < n_ijk[1])) || 
      //   ((pj[0] == 1) && (round_j[1] < n_ijk[1])));
      // testk = (round_k[1] >= k_idx[1]) &&  (((pk[1]<1) &&  (round_k[1] <= k_idx[le_k - 1])) || 
      //   ((pk[0]==1) && (round_k[0]<= k_idx[le_k - 1])));
      if (!testk & (pk[0]<1)) 
        testk = testk && (k_loc[round_k[1]] - k_loc[round_k[0]]== round_k[1]-round_k[0]);
      
      
      if (testi && testj && testk){
        
          
        
        //coin
        value[idx] = vol3D[(round_i[0] + (round_j[0] * n_ijk[0]) +  (k_loc[round_k[0]] * le_map))] * pi[0]*pj[0]*pk[0];
        if (pi[1]!=0) value[idx] = value[idx] + vol3D[(round_i[1] + (round_j[0] * n_ijk[0]) +  (k_loc[round_k[0]] * le_map))] * pi[1]*pj[0]*pk[0];
        if (pj[1]!=0) value[idx] = value[idx] + vol3D[(round_i[0] + (round_j[1] * n_ijk[0]) +  (k_loc[round_k[0]] * le_map))] * pi[0]*pj[1]*pk[0];
        if ((pi[1]!=0) && (pj[1]!=0)) value[idx] = value[idx] + vol3D[(round_i[1] + (round_j[1] * n_ijk[0]) +  (k_loc[round_k[0]] * le_map))] * pi[1]*pj[1]*pk[0];
        if (pk[1]!=0){
          value[idx] = value[idx] + vol3D[(round_i[0] + (round_j[0] * n_ijk[0]) +  (k_loc[round_k[1]] * le_map))] * pi[0]*pj[0]*pk[1];
          if (pi[1]!=0) value[idx] = value[idx] + vol3D[(round_i[1] + (round_j[0] * n_ijk[0]) +  (k_loc[round_k[1]] * le_map))] * pi[1]*pj[0]*pk[1]; 
          if (pj[1]!=0) value[idx] = value[idx] + vol3D[(round_i[0] + (round_j[1] * n_ijk[0]) +  (k_loc[round_k[1]] * le_map))] * pi[0]*pj[1]*pk[1];
          if ((pi[1]!=0) && (pj[1]!=0)) value[idx] = value[idx] + vol3D[(round_i[1] + (round_j[1] * n_ijk[0]) +  (k_loc[round_k[1]] * le_map))] * pi[1]*pj[1]*pk[1];
        }
        fill_i = (round_i[0] + 1) < round_i[1] ;
        fill_j = (round_j[0] + 1) < round_j[1] ;
        fill_k = (round_k[0] + 1) < round_k[1] ;
        //interieur
        
        if (fill_i && fill_j && fill_k){
          for(long i_ = round_i[0] + 1; i_ < round_i[1]; i_++)
            for( long j_ = round_j[0] + 1; j_< round_j[1]; j_++)
              for(long k_ = round_k[0] + 1; k_ < round_k[1]; k_++) 
                value[idx] = value[idx] + vol3D[(i_ + (j_ * n_ijk[0]) +  (k_loc[k_] * le_map))];
        }
        
        //arrete
        //i et j fixes 
        if (fill_k) for(long k_ = round_k[0] + 1; k_ < round_k[1]; k_++) {
          value[idx] = value[idx] + vol3D[(round_i[0] + (round_j[0] * n_ijk[0]) +  (k_loc[k_] * le_map))]*pi[0]*pj[0];
          if (pj[1]!=0) value[idx] = value[idx] + vol3D[(round_i[0] + (round_j[1] * n_ijk[0]) +  (k_loc[k_] * le_map))]*pi[0]*pj[1];
          if (pi[1]!=0) value[idx] = value[idx] + vol3D[(round_i[1] + (round_j[0] * n_ijk[0]) +  (k_loc[k_] * le_map))]*pi[1]*pj[0]; 
          if ((pi[1]!=0) && (pj[1]!=0)) value[idx] = value[idx] + vol3D[(round_i[1] + (round_j[1] * n_ijk[0]) +  (k_loc[k_] * le_map))]*pi[1]*pj[1]; 
        }
        
        //i et k fixes 
        if (fill_j) for( long j_ = round_j[0] + 1; j_< round_j[1]; j_++) {
          value[idx] = value[idx] +  vol3D[(round_i[0] + (j_ * n_ijk[0]) +  (k_loc[round_k[0]] * le_map))]*pi[0]*pk[0]; 
          if (pk[1]!=0) value[idx] = value[idx] + vol3D[(round_i[0] + (j_ * n_ijk[0]) +  (k_loc[round_k[1]] * le_map))]*pi[0]*pk[1]; 
          if (pi[1]!=0) value[idx] = value[idx] + vol3D[(round_i[1] + (j_ * n_ijk[0]) +  (k_loc[round_k[0]] * le_map))]*pi[1]*pk[0];
          if ((pi[1]!=0) && (pk[1]!=0)) value[idx] = value[idx] + vol3D[(round_i[1] + (j_ * n_ijk[0]) +  (k_loc[round_k[1]] * le_map))]*pi[1]*pk[1]; 
        }
        //j et k fixes 
        if (fill_i) for(long i_ = round_i[0] + 1; i_ < round_i[1]; i_++) {
          value[idx] = value[idx] + vol3D[(i_ + (round_j[0] * n_ijk[0]) +  (k_loc[round_k[0]] * le_map))]*pj[0]*pk[0];  
          if (pk[1]!=0) value[idx] = value[idx] + vol3D[(i_ + (round_j[0] * n_ijk[0]) +  (k_loc[round_k[1]] * le_map))]*pj[0]*pk[1];  
          if (pj[1]!=0) value[idx] = value[idx] + vol3D[(i_ + (round_j[1] * n_ijk[0]) +  (k_loc[round_k[0]] * le_map))]*pj[1]*pk[0];  
          if ((pj[1]!=0) && (pk[1]!=0)) value[idx] = value[idx] + vol3D[(i_ + (round_j[1] * n_ijk[0]) +  (k_loc[round_k[1]] * le_map))]*pj[1]*pk[1];
        }
        value[idx] = value[idx]/ (s_ijk[0]*s_ijk[1]*s_ijk[2]);
      }
    }
  }
  return (value);
}
