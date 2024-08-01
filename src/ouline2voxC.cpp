#include "esplib.h"

struct pt_st {
  double vol_index;
  double dY;
  double SR;
  int j;
  int maxj;
};


bool pt_sorter(pt_st const& lpt, pt_st const& rpt) {
  return (lpt.vol_index < rpt.vol_index);
}
// [[Rcpp::export(name = ".ouline2voxC")]] 
std::vector <double> ouline2voxC (std::vector <double> p,
                                  std::vector  <int> n_ijk,
                                  std::vector  <int> k_idx,
                                  std::vector  <int> k_loc,
                                  std::vector <double> O_ijk,
                                  std::vector <double> lambda_max,
                                  int ntab) {//, bool df
  
  std::vector <pt_st> ptl(ntab);

  std::vector <double> data(n_ijk[0]*n_ijk[1]*n_ijk[2]);
  std::vector <double> M1(3),M2(3), dM(2), roundM(3);
  int m,n,nb,ui,next_n, sens,nb_out,maxj; 
  double dl, dum;
  
  nb_out =0;
  for (m=0; m< (int) lambda_max.size(); m++){
    R_CheckUserInterrupt();
    roundM[2] = (int)( O_ijk[3*m + 2] + 0.5);
    if ((roundM[2]>=k_idx[0]) && (roundM[2]<= k_idx[n_ijk[2]-1])){
      if ((int) roundM[2]==k_idx[(int) roundM[2]]){
        nb = 0;
        data[nb] = 0 ;
        // Rcout<<"\n\n"<<std::to_string(m)<<"-- " <<std::to_string(nb_out);
        
        
        sens = signC(p[3*m]);
        // Rcout<<" |pi"<<std::to_string(p[3*m])<<" |sens"<<std::to_string(sens)<<"|n "<<std::to_string(n)<<"| ";
        if (sens !=0) {
          dum = 0.0;
          n = (int)(O_ijk[3*m]+0.5);
          if (sens==-1) { if (n_ijk[0]<n) n=n_ijk[0]; 
          } else{if (0>n) n=0;}
          while((n<=n_ijk[0]) &&  (n>=0) && (dum <= lambda_max[m])){
            n= n+sens;
            nb++;
            // data[nb]=  (int)((n - 0.5*sens - O_ijk[3*m])*1.0e6/ p[3*m])/1.0e6;
            data[nb]=  (n - 0.5*sens - O_ijk[3*m])/ p[3*m];
            dum = data[nb];
            // Rcout<<"i"<<std::to_string(data[nb])<< " ";
          }
        }
        
        sens = signC(p[3*m + 1]);
        // Rcout<<" |pj"<<std::to_string(p[3*m+1])<<" |sens" <<std::to_string(sens)<<"|n "<<std::to_string(n)<<"| ";
        if (sens !=0) {
          dum = 0.0;
          n = (int)(O_ijk[3*m + 1]+0.5);
          if (sens==-1) { if (n_ijk[1]<n) n=n_ijk[1]; 
          } else{if (0>n) n=0;}
          while((n<=n_ijk[1]) &&  (n>=0) && (dum<= lambda_max[m])){
            n= n+sens;
            nb++;
            // data[nb] = (int)((n - 0.5*sens - O_ijk[3*m + 1])*1.0e6/ p[3*m + 1])/1.0e6;
            data[nb] = (n - 0.5*sens - O_ijk[3*m + 1])/ p[3*m + 1];
            dum = data[nb];
            
            // Rcout<<"j"<<std::to_string(data[nb])<< " ";
          }
        }
        
        if (nb>1) std::sort(&(data[0]), &(data[0]) + nb);
        ui = 0;
        for(n=1; n<= nb; n++){
          if (data[ui]!= data[n]) {ui++;}
          data[ui] = data[n];
          if (data[n]>lambda_max[m]) break;
        }
        nb = ui;
        // Rcout<<"\n";
        // for(n=0; n<= nb; n++ ) Rcout<<std::to_string(data[n])<< " ";
        // Rcout<<"\n";
        n=0;
        
        while (n<nb) {
          // Rcout << std::to_string(n+1)<<" "<< std::to_string(data[n])<<"\n";
          next_n = n + 1;
          // Rcout <<"\n"<<  std::to_string(n)<<  "/" <<std::to_string(nb-1)<< " -  ";
          M1[0] = O_ijk[3*m] + data[n] * p[3*m];
          M1[1] = O_ijk[3*m + 1] + data[n] * p[3*m + 1];
          
          if (data[next_n]>lambda_max[m]) data[next_n]= lambda_max[m];
          M2[0] = O_ijk[3*m] + data[next_n] * p[3*m];
          M2[1] = O_ijk[3*m + 1] + data[next_n] * p[3*m + 1];
          // M2[2] = O_ijk[3*m + 2] + data[next_n] * p[3*m + 2];
          // Rcout << std::to_string(M1[0])<<"|"<< std::to_string(M1[1]) <<"  "<< std::to_string(M2[0])<<"|"<< std::to_string(M2[1]);
          
          // roundM[2] =(int) ((M1[2] + M2[2]+1)*0.5);
          roundM[0] = (int) ((M1[0] + M2[0]+1)*0.5);
          roundM[1] = (int) ((M1[1] + M2[1]+1)*0.5);
          
          if ((roundM[0]>=0) && (roundM[1]>=0)  && (roundM[0]<n_ijk[0]) && (roundM[1]<n_ijk[1]) ){
            dl = data[next_n] - data[n];
            
            dM[0] = dl * p[3*m];
            dM[1] = dl * p[3*m + 1];
            
            if (dM[1]!=0.0){
              ptl[nb_out].vol_index = roundM[0] + roundM[1]  * n_ijk[0] + (k_loc[roundM[2]] * n_ijk[0] *  n_ijk[1]);
              ptl[nb_out].dY = dM[1];
              ptl[nb_out].SR = (0.5 - M1[0]*0.5  -  M2[0]*0.5 + roundM[0]) * dM[1];
              ptl[nb_out].j = roundM[1];
              ptl[nb_out].maxj = n_ijk[0] - 1 + roundM[1]  * n_ijk[0] + (k_loc[roundM[2]] * n_ijk[0] *  n_ijk[1]);
              // Rcout <<" ind:" << std::to_string(ptl[nb_out].vol_index )<<" SR:" << std::to_string(ptl[nb_out].SR )<<" maxj:" <<std::to_string(ptl[nb_out].maxj );
              nb_out=nb_out + 1;
            }
          }  
          n = next_n;
        }
      }
    }
  }
  nb = nb_out;
  // Rcout<<"\n"<<std::to_string(round(((double) nb_out)*100.0/((double) ntab))/100.0)<<"  "<<std::to_string(nb_out)<<"/"<<std::to_string(ntab)<<" ";
  // for(n=0; n< nb; n++ ) Rcout<<"ind: "<<std::to_string(ptl[n].vol_index )<<", dY: "<<
  //   std::to_string(ptl[n].dY )<<", Sr: "<<std::to_string(ptl[n].SR )<<", maxj: "<<std::to_string(ptl[n].maxj )<< "\n";
  
  if (nb_out>1) std::sort(&(ptl[0]), &(ptl[0]) + nb, pt_sorter);
  
  // Rcout<<"\n\n";
  // for(n=0; n< nb; n++ ) Rcout<<"ind: "<<std::to_string(ptl[n].vol_index )<<", dY: "<<
  //   std::to_string(ptl[n].dY )<<", Sr: "<<std::to_string(ptl[n].SR )<<", maxj: "<<std::to_string(ptl[n].maxj )<< "\n";
  // 
  // 
  nb_out = 0;
  // Rcout<<"\n\n";
  for(n=1; n< nb; n++){
    if (ptl[nb_out].vol_index != ptl[n].vol_index) {
      nb_out++;
      ptl[nb_out].SR = ptl[n].SR;
      ptl[nb_out].dY = ptl[n].dY;
    } else{
      ptl[nb_out].SR += ptl[n].SR;
      ptl[nb_out].dY += ptl[n].dY;
    }
    ptl[nb_out].vol_index = ptl[n].vol_index;
    ptl[nb_out].j = ptl[n].j;
    ptl[nb_out].maxj = ptl[n].maxj;
  }
  nb_out++;
  // sens = signC(ptl[0].SR); 
  sens = 1; // attention , ici hypothèse que tous les contours tournent dans le sens des aiguilles d'une montre
  
  // if (df){
  //   std::vector <double> outdf(nb_out*3);
  //   if (nb_out==0) return(std::vector <double>());
  //   for(n=0;n<nb_out;n++){
  //     
  //     outdf[3*n] = ptl[n].vol_index;
  //     outdf[3*n + 1] = ptl[n].SR*sens;
  //     outdf[3*n + 2] = ptl[n].dY*sens;
  //   }
  //   
  //   return(outdf);
  // }
  
  for (n=0; n< (int) data.size();n++) data [n]=0.0;
  nb =0;
  // Rcout<<"\n"<<std::to_string(nb_out)<<"\n";
  // for(n=0; n< nb_out; n++ ) Rcout<<std::to_string(ptl[n].vol_index )<<" "<<
  //   std::to_string(ptl[n].dY )<<" "<<std::to_string(ptl[n].SR )<< " "<<
  //     std::to_string(ptl[n].maxj)<<"\n";
  // 
  
  while(nb<nb_out){
    R_CheckUserInterrupt();
    dum = data[ptl[nb].vol_index];
    data[ptl[nb].vol_index] = dum + sens * ptl[nb].SR;
    // Rcout<<"\n"<<std::to_string(nb+1)<<" val :"<<std::to_string(data[ptl[nb].vol_index] )<<"\n   "<<
    //   std::to_string(ptl[nb].vol_index)<<" "<<std::to_string(ptl[nb].SR)<<" "<<std::to_string(ptl[nb].dY);
    dum =  dum + sens * ptl[nb].dY;
    if (fabs(dum)<=1e-6) dum = 0.0;
    if (dum != 0.0){
      maxj = ptl[nb].maxj;
      if (nb!=nb_out-1) { if (ptl[nb].j == ptl[nb + 1].j)  maxj = ptl[nb+1].vol_index;}
      // Rcout<<"\n"<<std::to_string(nb) <<"  jusqu'à " << std::to_string(maxj)<<" remplit par "<<std::to_string(dum)<<" "<<std::to_string(ptl[nb].maxj);
      for(n = ptl[nb].vol_index + 1; n <= maxj; n++) data[n] =  dum;

    }
    nb++;

  }
  return(data);
}
