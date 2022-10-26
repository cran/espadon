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
// [[Rcpp::export(name = ".medianfilterC")]] 
std::vector <double> medianfilterC (
    std::vector <double> vol3D,
    std::vector <int> n_ijk,
    std::vector  <int> ball_i,
    std::vector  <int> ball_j,
    std::vector  <int> ball_k
){
  // std::vector  <float> max_slice_value,
  
  long idx, le, le_map,le_i, nb;
  long  i, j, k, rest;
  
  long jdx, new_jdx;
  int new_k, new_j, new_i;

  le = (int) vol3D.size ();
  le_i = (int) ball_i.size ();
  
  le_map = n_ijk[0] * n_ijk[1];
  std::vector <double> vect(le_i*le_i*le_i);
  std::vector <double> pr(le);

  for (idx= 0; idx < le; idx ++){
    k = (int) (idx/le_map);
    rest = idx - k*le_map;
    j = (int) (rest / n_ijk[0]);
    i = rest - j * n_ijk[0];
    
    nb = 0;
    // Rcout<<std::to_string(i)<<" "<<std::to_string(j)<<" "<<std::to_string(k)<<"\n";
    for (jdx = 0; jdx < le_i; jdx ++){
      new_k = k+ball_k[jdx];
      new_j = j+ball_j[jdx];
      new_i = i+ball_i[jdx];
      
      if ((new_k < n_ijk[2]) && (new_k >= 0) &&
          (new_j < n_ijk[1]) && (new_j >= 0) &&
          (new_i < n_ijk[0]) && (new_i >= 0)){
        new_jdx =  new_i + (new_j * n_ijk[0]) +  (new_k * le_map);
        vect[nb] = vol3D[new_jdx];
        nb = nb + 1;
      } 
    }
    if (nb>1) std::sort(&(vect[0]), &(vect[0]) + nb);
     
    if (nb % 2 != 0) {
      pr[idx] = vect[(int)nb / 2];
    } else {
      pr[idx]  = (vect[(int) ((nb - 1) / 2)] + vect[nb / 2]) / 2.0;
    } 

  }
  return(pr);
}
