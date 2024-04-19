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
  

// [[Rcpp::export(name = ".ptnearestC")]] 
Rcpp::List ptnearestC (
  std::vector <double> pt1_x,
  std::vector <double> pt1_y, 
  std::vector <double> pt1_z,
  std::vector <double> pt2_x,
  std::vector <double> pt2_y, 
  std::vector <double> pt2_z,
  bool full_info){
  
  unsigned int pt1_le, i,j;
  pt1_le = (unsigned int) pt1_x.size ();
  std::vector <int> index(pt1_le),index_min(2);
  std::vector <double> d(pt1_le);
  double d0,dmin; 
  
  dmin = sqrt(pow(pt2_x[0]-pt1_x[0],2) + pow(pt2_y[0]-pt1_y[0],2) + pow(pt2_z[0]-pt1_z[0],2));
  index_min[0] = 1;index_min[1] = 1;
  
  for(i = 0; i < pt1_le; i++){
    R_CheckUserInterrupt();
    d[i]= sqrt(pow(pt2_x[0]-pt1_x[i],2) + pow(pt2_y[0]-pt1_y[i],2) + pow(pt2_z[0]-pt1_z[i],2));
    index[i] = 1;

    for(j = 1; j < (unsigned int) pt2_x.size (); j++){
      d0 = sqrt(pow(pt2_x[j]-pt1_x[i],2) + pow(pt2_y[j]-pt1_y[i],2) + pow(pt2_z[j]-pt1_z[i],2));
      if (d0 < d[i]) {d[i] = d0; index[i]= j + 1;}
      if (d0 < dmin) {dmin = d0; index_min[0] = i + 1; index_min[1] = j + 1;}
    }
  }
  if (full_info){
    Rcpp::List res = Rcpp::List::create(
      Rcpp::Named("pt.index") =index_min,
      Rcpp::Named("d.min") =dmin,
      Rcpp::Named("d") = d,
      Rcpp::Named("index2") = index);
    return(res);
  } else{
    Rcpp::List res = Rcpp::List::create(
      Rcpp::Named("pt.index") =index_min,
      Rcpp::Named("d.min") =dmin);
    return(res);
  }
  
  
  }
  
  
