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
#define pi  3.141592653589793238462
// [[Rcpp::export(name = ".fansphereC")]] 
std::vector <double> fansphereC (double angle) {
  
  int N, M_theta, m, Ncount, M_phi, n;
  double a, d_theta, d_phi, theta, phi;
    
  a = 4*pi/pow(angle*pi/180.0,2);
  N = (int) a;
  if (a-N>0) N=N+1;
  a = 4*pi/N;
  M_theta = round (pi/sqrt(a));
  d_theta = pi/M_theta;
  d_phi = a/d_theta;
  N = 0;  
  for (m=0; m<M_theta; m++)
    N = N + round(2*pi*sin( pi*(m + 0.5)/M_theta)/d_phi);
  std::vector <double> p(5*N);
  Ncount = 0;
    
   
   for (m=0; m<M_theta; m++){
    theta = pi*(m + 0.5)/M_theta;
    M_phi = round(2*pi*sin(theta)/d_phi);
    for (n=0; n<M_phi; n++){
      phi = 2*pi*n/M_phi;
      p[5*Ncount]=sin(theta)*cos(phi);
      p[5*Ncount + 1]=sin(theta)*sin(phi);
      p[5*Ncount + 2]=cos(theta);
      p[5*Ncount + 3]=theta;
      p[5*Ncount + 4]=phi;
      if(fabs(p[5*Ncount])<1e-6) p[5*Ncount] = 0.0;
      if(fabs(p[5*Ncount + 1])<1e-6) p[5*Ncount + 1] = 0.0;
      if(fabs(p[5*Ncount + 2])<1e-6) p[5*Ncount + 2] = 0.0;
      
      Ncount  = Ncount+1;
    }
  }
   if (Ncount==0) return(std::vector <double>());
   std::vector <double> pr = {p.begin() ,p.begin() + 5*Ncount};
   return(pr);
}


// //[[Rcpp::export(name = ".test")]]
// std::vector <double> test(std::vector <double> index){
//   std::vector <double> t(3);
// 
//   if (index[0] ==5) {t[0] = 1.0;
//   }else {t[0] = 2.0;}
// 
//   t[1] = index[1] + 2;
//   Rcout<<std::to_string(NA_REAL)<<std::to_string(index[1])<<"\n";
//   t[2] = R_NaN + 2.0;
//   return (t);
// }
