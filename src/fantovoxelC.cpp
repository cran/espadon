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
#define pi  3.141592653589793238462
// [[Rcpp::export(name = ".fantovoxelC")]] 
std::vector <double> fantovoxelC (std::vector <double> p,
                                       std::vector  <int> n_ijk,
                                       std::vector  <int> k_idx,
                                       std::vector  <int> k_loc,
                                       std::vector <double> O_ijk,
                                       std::vector <double> vol_data,
                                       bool att = false,
                                       bool vol_value_flag = false,
                                       double vol_value = 1.0) {

  int p_le = (int) (p.size ()/3);
  int Ncount;
  if (att){
    Ncount = (n_ijk[0] + n_ijk[1] + n_ijk[2] + 4) * p_le * 5;
  } else {
    Ncount = (n_ijk[0] + n_ijk[1] + n_ijk[2] + 4) * p_le * 4;
  }
  int k_idx_le = (int) k_idx.size ();
  int roundk,vol_index;
  std::vector <double> outdf(Ncount);
  std::vector <double> lambda(n_ijk[0] + n_ijk[1] + n_ijk[2] + 4);
  std::vector <double> M_ijk(3);
  int m,n,nb,next_n,nb_out; 
  double dum,cum, dl;

  nb_out =0;
  for (m=0; m< p_le; m++){
    cum = 0;
    lambda[0] = 0 ;
    nb=1;
    if (p[3*m]!=0) {
      for(n=0; n<=n_ijk[0]; n++){
        dum =  round((n - 0.5 - O_ijk[0])*1.0e6/ p[3*m])/1.0e6;
        if (dum>=0){lambda[nb] = dum; nb ++;}
      }
    }
    if (p[3*m + 1]!=0) {
      for(n=0; n<=n_ijk[1]; n++){
        dum =  round((n - 0.5 - O_ijk[1])*1.0e6/ p[3*m + 1])/1.0e6;
        if (dum>=0){lambda[nb] = dum;  nb ++;}
      }
    }
    if (p[3*m + 2]!=0) {
      for(n=k_idx[0]; n<=k_idx[k_idx_le - 1] + 1 ; n++){
        dum =  round((n - 0.5 - O_ijk[2])*1.0e6/ p[3*m + 2])/1.0e6;
        if (dum>=0){lambda[nb] = dum;  nb ++;}
      }
    }
    if (nb>1) std::sort(&(lambda[0]), &(lambda[0]) + nb);
    n=0;

    while (n<nb - 1) {
      // Rcout << std::to_string(n+1)<<" "<< std::to_string(lambda[n])<<"\n";
      next_n = n + 1;
      while(lambda[next_n]==lambda[n] && next_n<nb-1) next_n++;
      if (lambda[next_n] == lambda[n]) break;
      M_ijk[0] = O_ijk[0] + 0.5 * (lambda[n] + lambda[next_n]) * p[3*m];
      M_ijk[1] = O_ijk[1] + 0.5 * (lambda[n] + lambda[next_n]) * p[3*m + 1];
      M_ijk[2] = O_ijk[2] + 0.5 * (lambda[n] + lambda[next_n]) * p[3*m + 2];

      if ((M_ijk[0] >= -0.5) && (M_ijk[1] >= -0.5) && (M_ijk[2] >= k_idx[0]-0.5) &&
        (M_ijk[0] < n_ijk[0] - 0.5) && (M_ijk[1] < n_ijk[1] - 0.5) &&
        (M_ijk[2] < k_idx[k_idx_le - 1] + 0.5)){

        roundk =(int) (M_ijk[2] + 0.5);
        if (roundk==k_idx[roundk]){
          vol_index = (int) (M_ijk[0] + 0.5) + ((int) (M_ijk[1] + 0.5)) * n_ijk[0] +
            (k_loc[roundk] * n_ijk[0] *  n_ijk[1]);
          dl = lambda[next_n] - lambda[n];
          cum = cum + (vol_data[vol_index] * dl); 
          if ((!vol_value_flag) || (vol_data[vol_index] == vol_value)){
            outdf[nb_out] = m + 1;
            outdf[nb_out + 1] = (double) vol_index + 1.0;
            outdf[nb_out + 2] = lambda[n];
            outdf[nb_out + 3] = dl;
            if (att) {
              outdf[nb_out + 4] = cum;
              nb_out=nb_out + 5;
            } else {nb_out=nb_out + 4;}
            
          }
        }
      }

      n = next_n;
    }

  }
  // return(lambda);
  if (nb_out==0) return(std::vector <double>());
  std::vector <double> pr = {outdf.begin() ,outdf.begin() + nb_out};
  return(pr);
  
}
