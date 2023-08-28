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

// [[Rcpp::export(name = ".isoclineC")]] 
std::vector <int> isoclineC (
    std::vector <int> it1,
    std::vector <int> it2){
  
  int nb_tt, le_it, i;
  int l_idx, pt_nb;
  int row_idx0,row_idx;
  int next_idx, where, dum;
  
  le_it = it1.size ();
  // std::vector <int> dispo_it (le_it);
  // std::vector <int> dispo_pt1 (le_it);
  // std::vector <int> dispo_pt2 (le_it);
  std::vector <int> dispo (4*le_it);
  for (i=0;i<le_it;i++)  {
    // dispo_it[i] = 0;
    // dispo_pt1[i] = 0;
    // dispo_pt2[i] = 0;
    dispo[4*i]=0;
    dispo[4*i + 1]=0;
    dispo[4*i + 2]=0;
    dispo[4*i + 3]=0; //swap 
  }
  nb_tt = 0;
  l_idx = 0;
  while(nb_tt< le_it){
    l_idx ++;
    pt_nb = 1;
    row_idx =-1;
    row_idx0=-1;
    for (i=0;i<le_it;i++){
      // if ( dispo_it[i] ==0) {
      if ( dispo[4*i] ==0) {
        row_idx = i;
        row_idx0 = i;
        break;
        }
    }

    next_idx = -1; 
    where = -1;
    
    while (row_idx!=-1){
      dispo[4*row_idx] = l_idx;
      dispo[4*row_idx + 1] = pt_nb;
      
      // dispo_it[row_idx] = l_idx;
      // dispo_pt1[row_idx] = pt_nb;
      nb_tt ++;
      if ((row_idx == row_idx0) && pt_nb!=1) break;

      if (it2[row_idx]!=-1) {
        pt_nb ++;
        dispo[4*row_idx + 2] = pt_nb;
        // dispo_pt2[row_idx] = pt_nb;
        for (i=0;i<le_it;i++){
          if ((dispo[4*i] == 0) && (it1[i] == it2[row_idx])) {
          // if ((dispo_it[i] == 0) && (it1[i] == it2[row_idx])) {
            next_idx = i;
            where = 1;
            break;
          }
          if ((dispo[4*i] == 0) && (it2[i] == it2[row_idx])) {
          // if ((dispo_it[i] == 0) && (it2[i] == it2[row_idx])) {
            next_idx = i;
            where = 2;
            break;
          }
        }
        if (where==2) {
          dum = it1[next_idx];
          it1[next_idx] = it2[next_idx];
          it2[next_idx] = dum;
          dispo[4*next_idx + 3] = 1;
        }
      }
      row_idx = next_idx;
      next_idx = -1; 
      where = -1;
    }
  }
  return(dispo);
}
