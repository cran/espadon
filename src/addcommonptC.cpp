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


// [[Rcpp::export(name = ".addcommonptC")]]   
Rcpp::List addcommonptC (
    std::vector <double> pt1_x,
    std::vector <double> pt1_y, 
    std::vector <double> pt2_x,
    std::vector <double> pt2_y, 
    std::vector <double> u1_x,
    std::vector <double> u1_y, 
    std::vector <double> u2_x,
    std::vector <double> u2_y,
    std::vector <double> d1,
    std::vector <double> d2,
    double eps = 1.0e-9){
  
  unsigned short pt1_le, pt2_le, i,j;
  unsigned short out1_i, out2_i;
  unsigned short out_i;
  Rcpp::List res;
  
  pt1_le = (unsigned short) pt1_x.size ();
  pt2_le = (unsigned short) pt2_x.size ();
  double A1A2x, A1A2y, u1u2, diff, proj1, proj2, k1l, k2l;//,k1l_, k2l_;
  
  std::vector <unsigned short> seg1((pt2_le-1)*(pt1_le-1)), seg2((pt2_le-1)*(pt1_le-1));
  std::vector <double> k1((pt2_le-1)*(pt1_le-1)), k2((pt2_le-1)*(pt1_le-1));

  out_i=0;
  // out1_i=0; out2_i=0;
  for (i = 0; i < pt1_le - 1; i++){
    R_CheckUserInterrupt();
    for(j = 0; j < pt2_le - 1; j++){
     
      u1u2 = u1_x[i]*u2_x[j] + u1_y[i] * u2_y[j];
      diff = 1-u1u2*u1u2;
      A1A2x = pt2_x[j]-pt1_x[i];
      
      A1A2y = pt2_y[j]-pt1_y[i];
      // Rcout<<"\n"<<std::to_string(i+1)<< " "<<std::to_string(j+1) <<" ";
      if ((A1A2x!=0.0) || (A1A2y!=0.0)){ // pas le même point
        
        if (fabs(diff)<eps){ // parallel
          // if (diff==0){ //
          // Rcout<<"parallel ";
          
          // if (fabs(A1A2x * u1_y[i]-u1_x[i]* A1A2y)<eps) { // same line
          if (A1A2x * u1_y[i]-u1_x[i]* A1A2y==0) { // same line
            // Rcout<<"same ";
            //  A1
            k1l = A1A2x * u1_x[i] + A1A2y * u1_y[i]; 
            
            if (fabs(k1l)<eps) k1l = 0.0;
            if (fabs(d1[i]-k1l)<eps) k1l = d1[i];
            
            if ((k1l>=0.0) && (k1l<d1[i])){
              // Rcout<<"add1";
              k1[out_i] = k1l;
              seg1[out_i] = i + 1;
              k2[out_i] = 0;
              seg2[out_i] = j + 1;
              out_i++;
            }
            
            k2l = -A1A2x * u2_x[j] - A1A2y * u2_y[j];
            
            if (fabs(k2l)<eps) k2l = 0.0;
            if (fabs(d2[j]-k2l)<eps) k2l = d2[j];
            
            if ((k2l>=0.0) && (k2l<d2[j])){
              // Rcout<<"add2";
              k1[out_i] = 0;
              seg1[out_i] = i + 1;
              k2[out_i] = k2l;
              seg2[out_i] = j + 1;
              out_i++;
            }
          } else{
            k1l = d1[i] + 1;
            k2l = d2[j] + 1;
          }
        } else{//secant
          

          proj1 = A1A2x * u1_x[i] + A1A2y * u1_y[i] ;
          proj2 = A1A2x * u2_x[j] + A1A2y * u2_y[j] ;
          k1l = (proj1 - proj2 * u1u2)/diff;
          k2l = (proj1 * u1u2 - proj2)/diff;
          // Rcout<<"secant; "<<std::to_string(k1l)<<"/"<<std::to_string(d1[i])<<";"<<
          //   std::to_string(k2l)<<"/"<<std::to_string(d2[j]);
          if (fabs(k1l)<eps) k1l = 0.0;
          if (fabs(k2l)<eps) k2l = 0.0;
          if (fabs(d1[i]-k1l)<eps) k1l = d1[i];
          if (fabs(d2[j]-k2l)<eps) k2l = d2[j];
          
          if ( (k1l>=0.0) && (k1l<d1[i])  && (k2l>=0.0) && (k2l<d2[j]) ){
            // Rcout<<"add";
            k1[out_i] = k1l;
            seg1[out_i] = i + 1;
            k2[out_i] = k2l;
            seg2[out_i] = j + 1;
            out_i++;
          }
          
        }
      }
    }
  }
  
  out1_i = out_i;
  out2_i = out_i;
  if ((out1_i==0)  && (out2_i==0)) {return(res);
  } else if ((out1_i!=0)  && (out2_i==0)){
    return(Rcpp::List::create(
        Rcpp::Named("pt1.index") = std::vector <unsigned short> {seg1.begin() ,seg1.begin() + out1_i},
        Rcpp::Named("k1") = std::vector <double> {k1.begin() ,k1.begin() + out1_i}));
  } else if ((out1_i==0)  && (out2_i!=0)){
    return(Rcpp::List::create(
        Rcpp::Named("pt2.index") = std::vector <unsigned short> {seg2.begin(), seg2.begin() + out2_i},
        Rcpp::Named("k2") = std::vector <double> {k2.begin() ,k2.begin() + out2_i}));
    
  } else {
    return(Rcpp::List::create(
        Rcpp::Named("pt1.index") = std::vector <unsigned short> {seg1.begin() ,seg1.begin() + out1_i},
        Rcpp::Named("k1") = std::vector <double> {k1.begin() ,k1.begin() + out1_i},
        Rcpp::Named("pt2.index") = std::vector <unsigned short> {seg2.begin(), seg2.begin() + out2_i},
        Rcpp::Named("k2") = std::vector <double> {k2.begin() ,k2.begin() + out2_i}));
  }
  
  
}

