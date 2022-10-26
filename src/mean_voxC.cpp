#include <Rcpp.h>
using namespace Rcpp;

#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>

#include <stdio.h>
#include <string.h>

#include <locale> 
// [[Rcpp::export(name = ".mean_voxC")]] 
std::vector <double> mean_voxC (
    std::vector <double> vol3D,
    std::vector <int> index,
    std::vector <int> index_list,
    std::vector <double> value_list,
    std::vector <double> value_att_list){
  
  int i;
  int le = vol3D.size ();
  
  std::vector <double> vol3D_som(le); 
  std::vector <double> vol3D_nb(le); 
  for (i=0;i<le;i++){
    vol3D_som[i] = vol3D[i];
    vol3D_nb[i] = 0;
  }
  for (i=0;i<(int) index_list.size ();i++) {
    vol3D_som[index_list[i]] = vol3D_som[index_list[i]] + value_list[i]*value_att_list[i];
    vol3D_nb[index_list[i]] =  vol3D_nb[index_list[i]]+ 1;
  }
  for (i=0;i<(int) index.size ();i++) {
    vol3D_som[index[i]] = vol3D_som[index[i]] / vol3D_nb[index[i]];
  }
  return(vol3D_som);
}
