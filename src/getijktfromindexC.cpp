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

// [[Rcpp::export(name = ".getijktfromindexC")]]
std::vector <int> getijktfromindexC(std::vector <int> index,
                                    std::vector  <int> k_idx,
                                    std::vector  <int> n_ijk){
  
  int le = (int) index.size ();
  int le_map = n_ijk[0] * n_ijk[1];
  int le_vol = le_map * n_ijk[2];
  int rest;
  
  std::vector <int> ijk(4*le);
  int i,Cindex;
  for(i=0; i<le; i++){
    Cindex = index[i] - 1;
    if ((Cindex<0)  || (Cindex>=le_vol)){
      ijk[4*i] =NA_INTEGER;
      ijk[4*i + 1] = NA_INTEGER;
      ijk[4*i + 2] = NA_INTEGER;
      ijk[4*i + 3] = 1;
    } else {
      ijk[4*i + 2] = (int) (Cindex / le_map);
      rest = (Cindex - ijk[4*i + 2] * le_map); 
      ijk[4*i + 1] = (int) (rest / n_ijk[0]);
      ijk[4*i] = rest - ijk[4*i + 1] *  n_ijk[0];
      ijk[4*i + 2] = k_idx[ijk[4*i + 2]];
      ijk[4*i + 3] = 1;
    }
    
  }
  return(ijk);
}
