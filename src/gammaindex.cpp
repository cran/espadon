#include <Rcpp.h>
using namespace Rcpp;

#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>
// #include <iostream>
#include <stdio.h>
#include <string.h>
#include <math.h>
  //#include <cstdint>
#include <locale> 

  // [[Rcpp::export(name = ".gammaindex")]] 
  std::vector <double> gammaindex (
      std::vector <double> vol3D,
      std::vector <double> vol3D_ref,
      std::vector <int> inspect_idx,
      std::vector <int> n_ijk,
      std::vector <double> rel_dxyz,
      std::vector  <int> ball_i,
      std::vector  <int> ball_j,
      std::vector  <int> ball_k,
      int around_idx,
      std::vector  <double> distance,
      double D_norm,
      bool local,
      double local_th_pc,
      double ref_pc ){
    
    bool find_gamma;
    long index,le, le_map, volindex;
    long  i, j, k, rest;
    int jdx, new_volindex, new_k, new_j, new_i;
    int jdx2, new_volindex2, new_k2, new_j2, new_i2, d2;
    double signD, curr_diff, diff_D1,diff_D2, th;
    double slope, dum;
    double curr_distance;
    double local_th = local_th_pc * D_norm;
    double a,b,dmin;
    
    le = (int) inspect_idx.size ();
    le_map = n_ijk[0] * n_ijk[1];
    std::vector <double> gamma(le);
    
    for (index= 0; index < le; index ++){
      // if (vol3D[idx]<=suppress_value_max)  {
      //   gamma[idx] = 0;
      // } else {
      R_CheckUserInterrupt();

      volindex = inspect_idx[index];
      k = (int) (volindex/le_map);
      rest = volindex - k*le_map;
      j = (int) (rest / n_ijk[0]);
      i = rest - j * n_ijk[0];

      if (local) {
        if (vol3D_ref[volindex]>=local_th) {th = vol3D_ref[volindex] * ref_pc;} 
        else {th =   ref_pc * local_th;}
      } else {th = ref_pc * D_norm;}
      
      // Rcout<<std::to_string(i)<<" "<<std::to_string(j)<<" "<<std::to_string(k)<<"\n";
      curr_diff = (vol3D[volindex] - vol3D_ref[volindex]);
     
      
      gamma[index] = fabs(curr_diff/th);
      find_gamma = false;
      if (curr_diff!= 0) {
        if (curr_diff>0) {signD = -1;} else {signD = 1;}
        for (jdx = 0; jdx < (int) ball_i.size (); jdx ++){
          new_k = k+ball_k[jdx]; new_j = j+ball_j[jdx]; new_i = i+ball_i[jdx];

          if ((new_k < n_ijk[2]) && (new_k >= 0) && (new_j < n_ijk[1]) &&
              (new_j >= 0) && (new_i < n_ijk[0]) && (new_i >= 0)){
            new_volindex =  new_i + (new_j * n_ijk[0]) +  (new_k * le_map);
            if (find_gamma && distance [jdx]> curr_distance) break;

            diff_D1 = vol3D[new_volindex] - vol3D_ref[volindex];
            if (signD*diff_D1>=0) {
              // on interpole :
              for (jdx2 = 0; jdx2 < around_idx; jdx2++){//on regarde su les 6 faces
                new_k2 = new_k+ball_k[jdx]; 
                new_j2 = new_j+ball_j[jdx]; 
                new_i2 = new_i+ball_i[jdx];

                if ((new_k2 < n_ijk[2]) && (new_k2 >= 0) && (new_j2 < n_ijk[1]) &&
                    (new_j2 >= 0) && (new_i2 < n_ijk[0]) && (new_i2 >= 0)){

                  new_volindex2 =  new_i2 + (new_j2 * n_ijk[0]) +  (new_k2 * le_map); 
                  d2 = sqrt(pow(rel_dxyz[0] * new_i2,2) + pow(rel_dxyz[1] * new_j2,2) +
                    pow(rel_dxyz[2] * new_k2,2)); // distance au point ref
                  diff_D2 = vol3D[new_volindex2] - vol3D_ref[volindex];
                  if ((signD*diff_D2<0) &&  (d2<distance[jdx])) {
                    slope = ((diff_D2-diff_D1)/th)/(d2 -distance[jdx]);
                    a = (1+pow(slope,2));
                    b = (diff_D1/th)- slope*distance[jdx];  
                    dmin = -b*slope/a;
                    dum =sqrt(pow(b,2)/(1+pow(slope,2)));
                    if(dum < gamma[index] && dmin<=distance[jdx] && dmin> d2) {
                      gamma[index] =dum; curr_distance = distance[jdx]; find_gamma = true;
                      }
                  }
                }
              }
             
              // Rcout<<std::to_string(gamma[index])<<" ";

            }else if (!find_gamma){
              dum = sqrt(pow(diff_D1 /th,2) + pow(distance[jdx],2));
              if (dum <=gamma[index]){
                gamma[index]  = dum;
                curr_distance = distance[jdx];
              } else if (distance [jdx] > gamma[index]) break;
            }
          }
        }
      }
    }
      
    return (gamma);
  }
  
    
  
