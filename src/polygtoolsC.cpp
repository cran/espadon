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







std::vector <double> polygsortC (std::vector <double> pol,
                                 int ncol = 3, 
                                 int compare =2, 
                                 bool clockwise= true){
  int le= (int)pol.size();
  int leord=(int)(le/ncol);
  if ((double)leord!=le/ncol) return(std::vector <double>());
  if (leord ==1) return(pol);
  
  int i, j, k;
  bool flag = true;
  for (i=0;i<ncol;i++) if (pol[i]!=pol[le-ncol+i]) {flag = false; break;}
  if (!flag) return(pol);
  int first=0;
  double isclw_v = (pol[0] - pol[(leord-2)*ncol] ) * (pol[ncol*1 + 1] - pol[1]) -   
    (pol[ncol*1] -pol[0]) * (pol[1] - pol[(leord-2)*ncol + 1]);
  double dum;
  std::vector <double> minp = {pol.begin() ,pol.begin() + compare};
  // Rcout<< "\n" << std::to_string(isclw_v)<< " ";
  for(i =1; i<leord-1; i++){
    dum = (pol[i*ncol] - pol[(i-1)*ncol] ) * (pol[ncol*(i+1) + 1] - pol[i*ncol + 1]) -   
      (pol[ncol*(i + 1)] -pol[i*ncol]) * (pol[i*ncol + 1] - pol[(i-1)*ncol + 1]);
    // Rcout<<std::to_string(dum)<<" ";
    for(j=0; j< compare; j++){
      
      if (pol[ncol*i+j]<minp[j]) {
        isclw_v = dum;
        for(k = j; j< compare; j++) minp[k] = pol[ncol*i+k];
        first = i;
        break;
      } else {
        if (isclw_v==0) isclw_v=dum;
        if (pol[ncol*i+j]>minp[j])break;
        }
      
    }
  }
  flag = (isclw_v < 0) == clockwise;
  // Rcout<<std::to_string(isclw_v);
  if ((first==0) && (flag)) return(pol);
  std::vector <double> outpol(le);
  if (flag){
    for (i=first;i<leord; i++) for(j=0; j< ncol; j++) outpol[(i-first)*ncol + j] = pol[i*ncol + j];
    for (i=0;i<first;i++) for(j=0; j< ncol; j++) outpol[(leord +i-first)*ncol + j] = pol[(i+1)*ncol + j];
  } else{
    for (i=first;i>=0; i--) for(j=0; j< ncol; j++) outpol[(first-i)*ncol + j] = pol[i*ncol + j];
    for (i=leord-1;i>first;i--) for(j=0; j< ncol; j++) outpol[(leord +first-i)*ncol + j] = pol[(i-1)*ncol + j];
  }
  return(outpol);
}

// [[Rcpp::export(name = ".polygcleanC")]] 
std::vector <double> polygcleanC (std::vector <double> pol, int ncol = 3, int compare =2, bool clockwise= true,
                                  double eps = 1.0e-6, bool sort = true){
  int le= (int)pol.size();
  if (le < ncol) return(std::vector <double>());
  int i,j;
  bool differ;
  std::vector <bool> diff(le-ncol);
  std::vector <double> outpol(le);
  for (i = 0; i < le; i++) pol[i] = eps * round(pol[i]/ eps);
  
  if (le ==ncol) return(pol);
  for (i = 0; i < le-ncol; i++) diff[i] = ((pol[i + ncol] - pol[i] )!= 0.0);
  i = 0;
  int outi;
  for (outi=0;outi<ncol; outi++) outpol[outi] = pol[outi];
  
  while (i<le-ncol){
    differ = false;
    for (j=i;j<i+compare; j++) if (diff[j]){differ = true; break;}
    if (differ) {
      for (j=0;j<ncol; j++) outpol[outi+j]= pol[i+j+ncol];
      outi = outi + ncol;
    }
    i = i + ncol;
  }
  std::vector <double> polout = {outpol.begin() ,outpol.begin() + outi};
  if (!sort) return(polout);
  return(polygsortC (polout,ncol,compare,clockwise));

}
  
