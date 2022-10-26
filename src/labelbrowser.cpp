#include <Rcpp.h>

using namespace Rcpp;

//------------------------------------------------------------------------------
// [[Rcpp::export(name = ".labelbrowser")]]
std::vector <unsigned int> labelbrowser (
    std::vector <bool> vol3D,
    std::vector  <unsigned int> n_ijk
){
  unsigned int i_ind, j_ind, k_ind;
  unsigned int i;
  
  unsigned int idx;
  unsigned int le_vol3D, le_map;
  unsigned int under_map, under_line, under_pt;
  unsigned int min_label, min_idx, max_idx;
  unsigned int label_under_pt, label_under_line, label_under_map;
  
  le_vol3D = vol3D.size ();
  le_map = n_ijk[0] *  n_ijk[1];
  std::vector <unsigned int> label(le_vol3D);
  std::vector <unsigned int> label_last_pt(le_vol3D);
  
  for (i = 0; i < le_vol3D; i++) {
    R_CheckUserInterrupt();
    if (vol3D[i] == false) {
      label[i] = le_vol3D;
      label_last_pt [i] = le_vol3D;
    } else {
      label[i] = i;
      label_last_pt[i] = i;
      k_ind = (int) (i / (le_map));
      j_ind = (int) ((i - (le_map * k_ind))/ n_ijk[0]);
      i_ind = i - (le_map * k_ind) - (n_ijk[0] * j_ind);
      
      if (i_ind==0) {
        under_pt = le_vol3D;
      } else {
        under_pt = i_ind - 1  + (n_ijk[0] * j_ind) + (le_map * k_ind);
        if (!vol3D[under_pt]) under_pt = le_vol3D;
      }
      
      if (j_ind==0) {
        under_line = le_vol3D;
      } else {
        under_line = i_ind + (n_ijk[0] * (j_ind - 1)) + (le_map * k_ind);
        if (!vol3D[under_line]) under_line = le_vol3D;
      }
        
      if (k_ind==0) {
        under_map = le_vol3D;
      } else {
        under_map = i_ind + (n_ijk[0] * j_ind) + (le_map * (k_ind-1));
        if (!vol3D[under_map]) under_map = le_vol3D;
      }   
        
       // 
       // Rcout <<"\n" << " " << std::to_string (i) << " " << std::to_string (i_ind) << " " << std::to_string (j_ind) << " " <<
       //   std::to_string (k_ind) << " " <<
       //   std::to_string (under_pt) << " " << std::to_string (under_line) << " " << std::to_string (under_map);
      
      if (under_pt != le_vol3D) label_under_pt = label[under_pt];
      else label_under_pt = le_vol3D;

      if (under_line != le_vol3D) label_under_line = label[under_line];
      else label_under_line = le_vol3D;

      if (under_map != le_vol3D) label_under_map = label[under_map];
      else label_under_map = le_vol3D;
    
    
      min_label = le_vol3D;
      min_idx = le_vol3D;
      max_idx = 0;
    
      if (label_under_pt != le_vol3D)  min_label = label_under_pt;
    
      if (label_under_line != le_vol3D) {
    
        if (label_under_line < min_label) {
          min_label = label_under_line;
          if (label_under_pt < le_vol3D) {
            min_idx = label_under_pt;
            max_idx = label_last_pt[label_under_pt];
          }
        } else if (label_under_line > min_label){
          min_idx =label_under_line;
          max_idx = label_last_pt[label_under_line];
        }
      }
    
      if (label_under_map != le_vol3D) {
    
        if (label_under_map <  min_label) {
          min_label = label_under_map;
    
          if ((label_under_pt != le_vol3D) && (label_under_pt != min_label)) {
            min_idx = label_under_pt;
            max_idx = label_last_pt[label_under_pt];
          }
    
          if ((label_under_line != le_vol3D) && (label_under_line != min_label)) {
            if (label_under_line < min_idx) min_idx = label_under_line;
            if (label_last_pt[label_under_line] > max_idx) max_idx = label_last_pt[label_under_line];
          }
        } else if (label_under_map >  min_label){
          if (label_under_map < min_idx) min_idx = label_under_map;
          if (label_last_pt[label_under_map] > max_idx) max_idx = label_last_pt[label_under_map];
        }
    
      }
    
    
    
      if (min_label != le_vol3D) {
        if (min_idx!=le_vol3D){
          for (idx = min_idx; idx<=max_idx; idx++) {
            if (label[idx]!=le_vol3D) {
              if (label[idx] == label_under_pt) label[idx] = min_label;
              if (label[idx] == label_under_line) label[idx] = min_label;
              if (label[idx] == label_under_map) label[idx] = min_label;
            }
          }
        }
    
        label[i] = min_label;
        label_last_pt[min_label] = i;
      }
    
    }
  }
  return label;
}


