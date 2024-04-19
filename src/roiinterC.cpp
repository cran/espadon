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

struct inter_pt_st {
  double k1;
  double vin;
  double vout;
};

bool pt_sorter(inter_pt_st const& lpt, inter_pt_st const& rpt) {
  return (lpt.k1 < rpt.k1);
}
  // [[Rcpp::export(name = ".roiinterC")]]   
  std::vector <double> roiinterC (
      std::vector <double> pt1_x,
      std::vector <double> pt1_y,
      double u1_x,
      double u1_y,
      double d1,
      std::vector <double> pt2_x,
      std::vector <double> pt2_y,
      double eps =1e-9){
  //
  unsigned short pt1_le, pt2_le;
  unsigned short i,j;
  unsigned short pt_j,limit_i, max_col=0;

  pt1_le = (unsigned short) pt1_x.size ();
  pt2_le = (unsigned short) pt2_x.size ();
  double A1A2x, A1A2y, u1u2, diff, proj1, proj2, k1l, k2l,vpin;

  std::vector <double> limits((pt2_le-1)*pt1_le);
  double vp_in, vp_out;
  std::vector <inter_pt_st> inter_pt(pt2_le-1);
  bool inlimits;
  double u2_x,u2_y,d2;


  for (i = 0; i < pt1_le ; i++){
    R_CheckUserInterrupt();
    pt_j=0;
    vp_in = (pt2_x[pt2_le - 1]-pt2_x[pt2_le - 2])* u1_y - u1_x * (pt2_y[pt2_le - 1]-pt2_y[pt2_le - 2]);
    // Rcout<<"\n"<<std::to_string(i) << " "<<std::to_string(pt1_x[i])<< " "<<std::to_string(pt1_y[i]) <<"\n    ";
    for(j = 0; j < pt2_le - 1; j++){
      limits[(j*pt1_le) +i]= (pt1_x[i] * u1_x) + (pt1_y[i] * u1_y);

      u2_x = pt2_x[j+1]-pt2_x[j];
      u2_y = pt2_y[j+1]-pt2_y[j];
      d2 = sqrt((u2_x * u2_x) + (u2_y*u2_y));
      u2_x = u2_x/d2; u2_y = u2_y/d2;
      
      vp_out = u2_x * u1_y- u1_x * u2_y;
  
      // if ((pt2_x[j+1]-pt1_x[i]) &&  (pt2_y[j+1]!=pt1_y[i])) {
      if ((((pt2_x[j+1]-pt1_x[i])*(pt2_x[j]-pt1_x[i])<=0) && (u1_x==0.0)) ||
      (((pt2_y[j+1]-pt1_y[i])*(pt2_y[j]-pt1_y[i])<=0) && (u1_y==0.0))){
        
        u1u2 = u1_x*u2_x + u1_y * u2_y;
        diff = 1-u1u2*u1u2;
        A1A2x = pt2_x[j]-pt1_x[i];
        A1A2y = pt2_y[j]-pt1_y[i];
        // if ((A1A2x!=0.0) || (A1A2y!=0.0)){ // pas le même point
        if (fabs(diff)<=eps){ // parallel
          // if (diff==0){ // parallel
          if (fabs(A1A2x * u1_y-u1_x* A1A2y)<=eps) { // same line
            // if (A1A2x * u1_y-u1_x* A1A2y==0) { // same line
            k1l = A1A2x * u1_x + A1A2y * u1_y;
            // if (fabs(k1l)<=eps) k1l = 0.0;
            // if (fabs(d1-k1l)<=eps) k1l = d1;
            // Rcout<<"    "<<std::to_string((pt2_x[j+1]-pt1_x[i]))<<" "<<std::to_string((pt2_x[j]-pt1_x[i]))<<
            //   " "<<std::to_string((pt2_y[j+1]-pt1_y[i])) <<" "<<std::to_string(pt2_y[j]-pt1_y[i]);
            // if ((k1l>=0.0) && (k1l<d1)){
              inter_pt[pt_j].k1 = k1l;
              // Rcout<<"par("<<std::to_string(j)<<" "<<std::to_string(k1l + (pt1_x[i] * u1_x) + (pt1_y[i] * u1_y))<<" "<<std::to_string(u1u2)<<
              //   " "<<std::to_string(vp_in)<<" "<<std::to_string(vp_out)<<") ";
              if (u1u2>=0){
                inter_pt[pt_j].vin = vp_in;
                inter_pt[pt_j].vout = vp_out;
              } else{
                inter_pt[pt_j].vin = vp_out;
                inter_pt[pt_j].vout = vp_in;
              }
              pt_j++;
            // }
          }
        } else{//secant
          proj1 = A1A2x * u1_x + A1A2y * u1_y;
          proj2 = A1A2x * u2_x + A1A2y * u2_y;
          k1l = (proj1 - proj2 * u1u2)/diff;
          k2l = (proj1 * u1u2 - proj2)/diff;
          // if (fabs(k1l)<=eps) k1l = 0.0;
          if (fabs(k2l)<=eps) k2l = 0.0;
          // if (fabs(d1-k1l)<=eps) k1l = d1;
          if (fabs(d2-k2l)<=eps) k2l = d2;
          if ( (k1l>=0.0) && (k1l<d1)  && (k2l>=0.0) && (k2l<d2) ){
            // Rcout<<"    "<<std::to_string((pt2_x[j+1]-pt1_x[i]))<<" "<<std::to_string((pt2_x[j]-pt1_x[i]))<<
            //   " "<<std::to_string((pt2_y[j+1]-pt1_y[i])) <<" "<<std::to_string(pt2_y[j]-pt1_y[i]);
            inter_pt[pt_j].k1 = k1l;
            // Rcout<<"sec("<<std::to_string(j)<<" "<<std::to_string(k1l + (pt1_x[i] * u1_x) + (pt1_y[i] * u1_y))<<
            //   " "<<std::to_string(k2l)<< "|"<<std::to_string((d2[j]))<<" "<<std::to_string(u1u2)<<
            //     " "<<std::to_string(vp_in)<<" "<<std::to_string(vp_out)<<") ";
            if (k2l==0.0){
              if (u1u2>=0){
                inter_pt[pt_j].vin = vp_in;
                inter_pt[pt_j].vout = vp_out;
              } else{
                inter_pt[pt_j].vin = vp_out;
                inter_pt[pt_j].vout = vp_in;
              }
            }else {
              inter_pt[pt_j].vin = vp_out;
              inter_pt[pt_j].vout = vp_out;
            }
            pt_j++;
          }
        }
      }
      vp_in = vp_out;

    }

    if (pt_j>0){
      if (pt_j>1) std::sort(&(inter_pt[0]), &(inter_pt[0]) + pt_j, pt_sorter);



    limit_i =0;
    inlimits = false;

    if (pt_j>1){

      for (j=0; j<pt_j; j++) {
        if (inlimits){
          if (inter_pt[j].vout * vpin < 0) {
            inlimits = false;
            vpin = 0;
            limits[(limit_i * pt1_le) + i] += inter_pt[j].k1;
            // Rcout<<"("<<std::to_string ((limit_i * pt1_le) + i + 1)<<" "<<std::to_string (limits[(limit_i * pt1_le) + i])<<") ";
            limit_i++;
          }
        } else{
          limits[(limit_i * pt1_le) + i] += inter_pt[j].k1;
          // Rcout<<"("<<std::to_string ((limit_i * pt1_le) + i + 1)<<" "<<std::to_string (limits[(limit_i * pt1_le) + i])<<
            // " "<<std::to_string (inter_pt[j].vout)<<" "<<std::to_string (inter_pt[j].vin)<<") ";
          limit_i++;
          if (inter_pt[j].vout * inter_pt[j].vin  < 0){
            limits[(limit_i * pt1_le) + i] += inter_pt[j].k1;
            // Rcout<<"("<<std::to_string ((limit_i * pt1_le) + i + 1)<<" "<<std::to_string (limits[(limit_i * pt1_le) + i])<<") ";
            limit_i++;
          } else{
            inlimits = true;
            vpin = inter_pt[j].vin;
          }

        }
      }

    } else{
      limits[(limit_i * pt1_le) + i] += inter_pt[0].k1;
      // Rcout<<"\n   ("<<std::to_string ((limit_i * pt1_le) + i + 1)<<" "<<std::to_string (limits[(limit_i * pt1_le) + i])<<") ";
      limit_i++;
      limits[(limit_i * pt1_le) + i] += inter_pt[0].k1;
      // Rcout<<"\n   ("<<std::to_string ((limit_i * pt1_le) + i + 1)<<" "<<std::to_string (limits[(limit_i * pt1_le) + i])<<") ";
      limit_i++;
    }

    if (limit_i>max_col) max_col=limit_i;

    // for(j = 0; j<max_col+1;j++){
    //   Rcout<<"\n    "<<std::to_string (j)<<"   ";
    //   for(pt_j=0;pt_j<4;pt_j++) Rcout<<"("<< std::to_string ((j*pt1_le) + pt_j)<<" "<<std::to_string (limits[(j*pt1_le) + pt_j])<<")";
    // }
    // Rcout<<"\n    colonnes:"<<std::to_string (max_col)<<"\n";

    }

  }

  
  // return(limits);
  return(std::vector <double> {limits.begin() ,limits.begin() + max_col*pt1_le});
 
  
}
    
    
