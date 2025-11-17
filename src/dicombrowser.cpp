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

//Bit_Order
#define High_Order_First  1
#define Low_Order_First 0
#define default_endianness Low_Order_First;

// VRness
#define explicit_VR 0
#define implicit_VR 1
#define default_VRness  explicit_VR
#define Undef_Length 0xFFFFFFFF


// [[Rcpp::export(name = ".dicombrowser")]]
Rcpp::List dicombrowser (std::vector <unsigned char> dicomrawdata,  
                          DataFrame tagdico, unsigned int nbTAG=0,
                          std::string stop_tag="", unsigned int stop_level=0,
                          bool full_info = false,  bool verbose = false){
  
  
  std::list<std::string> vTAG;
  std::list<std::string> vVR;
  std::list<std::string> vendian;
  std::list<int> vstart;
  std::list<int> vstop;
  std::list<std::string> vencapsload;
  std::list<int> vloadstart;
  std::list<int> vloadstop;
  std::list<int> vtagstart;
  
  
  std::string list_el;
  const std::string special_VR[] = {"OB", "OD", "OW", "OF", "OL","SQ", "UC", "UR", "UT", "UN"};
  const unsigned int length_special_VR[]  = {1, 8, 2, 4, 1,1,1, 1, 1, 1};
  const int encaps_max = 20;
  
  bool preamble,find;
  unsigned long long int my_cursor, counter_start, my_VL_start, my_VL_stop;
  int VRness, decoded_VRness, endianness, decoded_endianness;
  const std::string endian[]= {"little","big"};
  unsigned int tag_number = 0;
  unsigned long long int length_of_Group;
  unsigned int encapsulation_level,encaps_dum;
  unsigned int curr_encaps =0;
  std::string Group_TAG; 
  std::string tag="(0000,0000)";
  std::string tag_hhhh_hhxx="(0000,00xx)";
  std::string tag_hhhh_hhxh="(0000,00x0)";
  std::string tag_hhhh_xxxh="(0000,xxx0)";
  std::string tag_hhhh_xxxx="(0000,xxxx)";
  std::string tag_hhxx_hhhh="(00xx,0000)";
  uint16_t grouptag, elementtag;
  uint32_t my_VL=0;
  std::string std_my_VR, stdgrouptag1,stdgrouptag, stdelementtag, dcm;
  std::string grouptag_le="0000";
  std::string character_set = "0123456789ABCDEF";
  unsigned int length_data=0;
  unsigned int encapsulation_load[encaps_max];
  unsigned int item_nb[encaps_max];
  std::string tag_name[encaps_max];
  
  unsigned int offset_count = 0;
  
  unsigned int item_index;
  CharacterVector col_tag = tagdico["tag"];
  CharacterVector col_VR = tagdico["VR"];
  
  
  // Initialisation
  //----------------
  preamble = true;
  if (dicomrawdata.size () < 132) {
    preamble = false;
  } else {
    for (int i = 0; i < 128;i++) {
      if (dicomrawdata[i] != 0) preamble = false;
    }
    if ((dicomrawdata[128] != 68) || (dicomrawdata[129] != 73) || (dicomrawdata[130] != 67) ||
        (dicomrawdata[131] != 77))  preamble = false;
  }
  
  for (int i=0; i < encaps_max; i++) {encapsulation_load[i] = 0; item_nb[i] = 0; tag_name[i] = "";}
  
  
  
  if (preamble==TRUE) {
    my_cursor = 132;
    VRness = default_VRness;
  } else {
    my_cursor = 0;
    VRness = implicit_VR;
  }
  
  endianness = default_endianness;
  decoded_endianness = default_endianness;
  decoded_VRness = VRness;
  encapsulation_level=0;
  length_of_Group=0;           
  // parcours
  //---------
  unsigned int event = 0;
  while (  (my_cursor < dicomrawdata.size ()) && 
           ((tag_number < nbTAG) || (nbTAG == 0)) &&
           ((stop_tag == "") || (tag < stop_tag) || (curr_encaps != stop_level))) {
    
    if (event % 10000 == 0) R_CheckUserInterrupt();
    
    if ((encapsulation_level > 0) && (encapsulation_load[encapsulation_level] == 0)) {
      encaps_dum = encapsulation_level;
      for (int i = encaps_dum; i > 0; i--) {
        if (encapsulation_load [i] == 0) {
          item_nb[encapsulation_level] = 0;
          tag_name[encapsulation_level] = "";
          encapsulation_level = encapsulation_level - 1;
        }
        else break;
      }
    } else {
      
      item_nb[encapsulation_level] = item_nb[encapsulation_level] + 1;
      counter_start = my_cursor;
      
      // recuperation du TAG
      if (endianness == Low_Order_First) {
        grouptag = (uint16_t) dicomrawdata [my_cursor] + (256 * (uint16_t) dicomrawdata [my_cursor + 1]) ;
        elementtag  = (uint16_t) dicomrawdata [my_cursor + 2] + (256 * (uint16_t) dicomrawdata [my_cursor + 3]) ;
        
      } else {
        grouptag = (uint16_t) dicomrawdata [my_cursor + 1] + (256 * (uint16_t) dicomrawdata [my_cursor]) ;
        elementtag  = (uint16_t) dicomrawdata [my_cursor + 3] + (256 * (uint16_t) dicomrawdata [my_cursor + 2]) ;
      }
      
      stdgrouptag1 = {character_set.at((int)(grouptag / 0x1000)) ,character_set.at((int)((grouptag / 0x0100) % 0x10))};
      stdgrouptag = {character_set.at((int)((grouptag / 0x0010)  % 0x10)), character_set.at((int)(grouptag % 0x10))};
      stdelementtag = {character_set.at((int)(elementtag / 0x1000)) ,character_set.at((int)((elementtag / 0x0100) % 0x10)),
                       character_set.at((int)((elementtag / 0x0010)  % 0x10)), character_set.at((int)(elementtag % 0x10))};
      stdgrouptag =  stdgrouptag1 + stdgrouptag;
      tag = "(" + stdgrouptag + "," + stdelementtag + ")";
      tag_hhhh_hhxx = "(" + stdgrouptag + "," +  stdelementtag[0] +  stdelementtag[1]+ "xx)";
      tag_hhhh_hhxh = "(" + stdgrouptag + "," +  stdelementtag[0] +  stdelementtag[1]+ "x" +  stdelementtag[3] + ")";
      tag_hhhh_xxxh = "(" + stdgrouptag + ",xxx" +  stdelementtag[3] + ")";
      tag_hhhh_xxxx = "(" + stdgrouptag + "," +  "xxxx)";
      tag_hhxx_hhhh = "(" + stdgrouptag1 +  "xx," + stdelementtag + ")";
      
      curr_encaps = encapsulation_level;
      // placement aux octets suivants
      my_cursor  = my_cursor + 4;
      offset_count = 4;
      tag_number=tag_number+1;
      
      if ((tag_number<4) && (grouptag != 0x0002) && (grouptag != 0x0008)) {
        return(Rcpp::DataFrame::create(Rcpp::Named("error") = "not dicom compliant"));
      }
      
      
      if ((grouptag == 0xFFFE) && ((elementtag == 0xE0DD) || (elementtag == 0xE00D)) ){// Sequence_Delimitation_Tag ou Item_Delimitation_Tag
        //-----------------------------------------------------------------------------  
        tag_name[encapsulation_level] = tag;
        
        //        // affichage
        list_el = "";
        for (unsigned int i = 0; i< encapsulation_level; i++) {list_el = list_el+ tag_name[i] + " ";}
        vTAG.push_back(list_el + tag_name[encapsulation_level]);
        vVR.push_back("00");
        vendian.push_back(endian[endianness]);
        vstart.push_back(NA_INTEGER);
        vstop.push_back(NA_INTEGER);
        
        
        if (full_info){
          list_el = "";
          if (encapsulation_level>0) for (unsigned int i = 1; i< encapsulation_level+1; i++) {
            if (encapsulation_load[i] == Undef_Length) list_el = list_el +  "NA ";
            else list_el = list_el+ std::to_string (encapsulation_load[i]) + " ";
          }
          vencapsload.push_back(list_el);
          vtagstart.push_back(counter_start + 1);
          vloadstart.push_back(NA_INTEGER);
          vloadstop.push_back(NA_INTEGER);
        }
        // if (verbose)  Rcout << "\n" << list_el;
        
        
        
        my_cursor = my_cursor + 4;
        offset_count =offset_count + 4;
        for (unsigned int i = 0; i<=encapsulation_level; i++) {
          if ((encapsulation_load [i] != Undef_Length) && (encapsulation_load [i] != 0)) encapsulation_load [i] = encapsulation_load [i] - offset_count;
        }
        
        encapsulation_load [encapsulation_level] = 0;
        item_nb[encapsulation_level] = 0;
        tag_name[encapsulation_level]  = "";
        if (encapsulation_level>0)  encapsulation_level= encapsulation_level - 1;
        
        
      } else if ((grouptag == 0xFFFE) && (elementtag == 0xE000)){ // item tag
        //--------------------------------------------------------
        
        tag_name[encapsulation_level]= "item" + std::to_string (item_nb[encapsulation_level]);
        
        my_VL_start = my_cursor + 1;
        my_VL_stop = my_cursor + 4;
        if (endianness == Low_Order_First) my_VL = (uint32_t) dicomrawdata [my_cursor] + (256 * (uint32_t) dicomrawdata [my_cursor + 1]) +
          (65536 * (uint32_t) dicomrawdata [my_cursor + 2]) + (16777216 * (uint32_t) dicomrawdata [my_cursor +3]);
        else my_VL = (uint32_t) dicomrawdata [my_cursor + 3] +  (256 * (uint32_t) dicomrawdata [my_cursor + 2]) +
          (65536 * (uint32_t) dicomrawdata [my_cursor + 1]) + (16777216 * (uint32_t) dicomrawdata [my_cursor]);
        
        my_cursor = my_cursor + 4;
        offset_count = offset_count + 4;     
        
        //--
        list_el = "";
        for (unsigned int i = 0; i< encapsulation_level; i++) {list_el = list_el+ tag_name[i] + " ";}
        
        
        vTAG.push_back(list_el + tag_name[encapsulation_level]);
        vVR.push_back("00");
        vendian.push_back(endian[endianness]);
        if(my_VL == 0){
          vstart.push_back(NA_INTEGER);
          vstop.push_back(NA_INTEGER);
        } else if (my_VL == Undef_Length) {
          vstart.push_back(my_cursor + 1);
          vstop.push_back(NA_INTEGER);
        } else  {
          vstart.push_back(my_cursor + 1);
          vstop.push_back(my_cursor + my_VL);
        }
        
        
        
        if (full_info){
          list_el = "";
          if (encapsulation_level>0) for (unsigned int i = 1; i< encapsulation_level+1; i++) {
            if (encapsulation_load[i] == Undef_Length) list_el = list_el +  "NA ";
            else list_el = list_el+ std::to_string (encapsulation_load[i]) + " ";
          }
          vencapsload.push_back(list_el);
          vtagstart.push_back(counter_start + 1);
          vloadstart.push_back(my_VL_start);
          vloadstop.push_back(my_VL_stop);
        }
        
        //--            
        
        if ((my_VL != Undef_Length) && (my_cursor+my_VL>dicomrawdata.size())) {
          return(Rcpp::List::create(Rcpp::Named("error") = "not dicom compliant"));}
        
        for (unsigned int i = 0; i<=encapsulation_level; i++) {
          if ((encapsulation_load [i] != Undef_Length) && (encapsulation_load [i] != 0)) encapsulation_load [i] = encapsulation_load [i] - offset_count;
        }
        encapsulation_level = encapsulation_level + 1;
        encapsulation_load[encapsulation_level] = my_VL;
        
      } else { //all others cases
        //--------------------------
        tag_name[encapsulation_level] = tag;
        
        if (VRness == explicit_VR){
          std_my_VR = ""; std_my_VR.append (1, dicomrawdata[my_cursor]); std_my_VR.append (1, dicomrawdata[my_cursor+1]);
          my_cursor = my_cursor + 2;
          offset_count = offset_count + 2;  
        } else{
          if (stdgrouptag=="1010") {
            std_my_VR = "US";
          } else {
            find = false;
            for(item_index=0; item_index< (unsigned int)col_tag.length(); item_index++){
              
              if (col_tag[item_index]==tag) find=true;
              if (col_tag[item_index]==tag_hhhh_hhxx) find=true;
              if (col_tag[item_index]==tag_hhhh_hhxh) find=true;
              if (col_tag[item_index]==tag_hhhh_xxxh) find=true;
              if (col_tag[item_index]==tag_hhxx_hhhh) find=true;
              if (find) break;
            }
            if (!find) {
              if (tag_number < 4) {return(Rcpp::List::create(Rcpp::Named("error") = "not dicom compliant"));}
              if  (stdelementtag == "0000") {std_my_VR = "UL";
              } else {std_my_VR ="UN";}
            } else {
              std_my_VR = (std::string) CHARACTER_VALUE(col_VR[item_index]);
            }
          }
        }
        
        
        length_data = 0;
        for (int i=0; i < 10; i++) if (std_my_VR == special_VR[i]) {length_data = length_special_VR[i]; break;}
        
        if ((length_data !=0) && (VRness == explicit_VR)) {
          my_cursor = my_cursor + 2;
          offset_count =offset_count+2;
        }
        if (((length_data !=0) && (VRness == explicit_VR)) || (VRness != explicit_VR)){
          my_VL_start = my_cursor + 1;
          my_VL_stop = my_cursor + 4;
          
          if (endianness == Low_Order_First) my_VL = (uint32_t) dicomrawdata [my_cursor] + (256 * (uint32_t) dicomrawdata [my_cursor + 1]) +
            (65536 * (uint32_t) dicomrawdata [my_cursor + 2]) + (16777216 * (uint32_t) dicomrawdata [my_cursor +3]);
          else my_VL = (uint32_t) dicomrawdata [my_cursor + 3] +  (256 * (uint32_t) dicomrawdata [my_cursor + 2]) +
            (65536 * (uint32_t) dicomrawdata [my_cursor + 1]) + (16777216 * (uint32_t) dicomrawdata [my_cursor]);
          my_cursor = my_cursor + 4;
          offset_count = offset_count + 4;
        } else {
          my_VL_start = my_cursor + 1;
          my_VL_stop = my_cursor + 2;
          
          if (endianness == Low_Order_First)  my_VL = (uint32_t) dicomrawdata [my_cursor] + (256 * (uint32_t) dicomrawdata [my_cursor + 1]);
          else  my_VL = (uint32_t) dicomrawdata [my_cursor + 1] +  (256 * (uint32_t) dicomrawdata [my_cursor]);
          my_cursor = my_cursor + 2;
          offset_count = offset_count + 2;
        }
        
        // --
        list_el = "";
        for (unsigned int i = 0; i< encapsulation_level; i++) {list_el = list_el+ tag_name[i] + " ";}
        vTAG.push_back(list_el + tag_name[encapsulation_level]);
        vVR.push_back(std_my_VR);
        vendian.push_back(endian[endianness]);
        
        if(my_VL == 0){
          vstart.push_back(NA_INTEGER);
          vstop.push_back(NA_INTEGER);
        } else if (my_VL == Undef_Length) {
          vstart.push_back(my_cursor + 1);
          vstop.push_back(NA_INTEGER);
        } else  {
          vstart.push_back(my_cursor + 1);
          vstop.push_back(my_cursor + my_VL);
        }
        
        
        if (full_info){
          list_el = "";
          if (encapsulation_level>0) for (unsigned int i = 1; i< encapsulation_level+1; i++) {
            if (encapsulation_load[i] == Undef_Length) list_el = list_el +  "NA ";
            else list_el = list_el+ std::to_string (encapsulation_load[i]) + " ";
          } 
          
          vencapsload.push_back(list_el);
          vtagstart.push_back(counter_start + 1);
          vloadstart.push_back(my_VL_start);
          vloadstop.push_back(my_VL_stop);
        }
        
        // if (verbose)  Rcout << "\n" << list_el;
        // --
        
        if ((my_VL != Undef_Length) && (my_cursor+my_VL>dicomrawdata.size())) {
          return(Rcpp::List::create(Rcpp::Named("error") = "not dicom compliant"));}
        
        // decode endianness et VRness
        if ((stdgrouptag == "0002") && (stdelementtag == "0010")){
          dcm = "";
          if (my_VL>0) for (uint32_t i = 0; i < my_VL-1; i++) dcm.append (1, dicomrawdata[my_cursor + i]);
          if (dicomrawdata [my_cursor + my_VL-1] != 0)  dcm.append (1, dicomrawdata[my_cursor + my_VL-1]);
          
          if (dcm =="1.2.840.10008.1.2") {
            decoded_endianness = Low_Order_First;
            decoded_VRness = implicit_VR;
            if (verbose) Rcout << "low order - implicit\n";
          } else if  (dcm == "1.2.840.10008.1.2.1") {
            decoded_endianness = Low_Order_First;
            decoded_VRness = explicit_VR;
            if (verbose) Rcout << "low order - explicit\n";
          } else if  (dcm == "1.2.840.10008.1.2.1.99") {
            decoded_endianness = Low_Order_First; // deflated ???
            decoded_VRness = explicit_VR;
            if (verbose) Rcout << "low order - explicit\n";
          } else if  (dcm == "1.2.840.10008.1.2.2") {
            decoded_endianness = High_Order_First;
            decoded_VRness = explicit_VR;
            if (verbose) Rcout << "high order - explicit\n";
          }  
        }
        
        // fin decode endianness et VRness
        
        // recuperation de la longueur du header
        if ( (stdgrouptag == "0002") &&  (stdelementtag == "0000")) {
          // if  (stdelementtag == "0000") {
          grouptag_le = stdgrouptag; 
          if (endianness == Low_Order_First) length_of_Group = (uint32_t) dicomrawdata [my_cursor] + (256 * (uint32_t) dicomrawdata [my_cursor + 1]) +
            (65536 * (uint32_t) dicomrawdata [my_cursor + 2]) + (16777216 * (uint32_t) dicomrawdata [my_cursor +3]);
          else length_of_Group = (uint32_t) dicomrawdata [my_cursor + 3] +  (256 * (uint32_t) dicomrawdata [my_cursor + 2]) +
            (65536 * (uint32_t) dicomrawdata [my_cursor + 1]) + (16777216 * (uint32_t) dicomrawdata [my_cursor]);
        }
        // recuperation de la longueur du header
        
        if (std_my_VR == "AT") {
          my_cursor = my_cursor +  4;
          offset_count = offset_count + 4;
        } else if (!((std_my_VR=="SQ") || (std_my_VR == "00"))){
          my_cursor = my_cursor + my_VL;
          offset_count = offset_count + my_VL;
        }
        // prise en compte de endianness et VRness
        if ((stdgrouptag == "0002")  &&  (stdelementtag != "0000")){
          length_of_Group = length_of_Group - (my_cursor - counter_start);
          if (length_of_Group == 0) {
            VRness = decoded_VRness;
            endianness = decoded_endianness;
          }
        }
        // fin prise en compte de endianness et VRness
        
        for (unsigned int i = 0; i<=encapsulation_level; i++) {
          if ((encapsulation_load [i] != Undef_Length) && (encapsulation_load [i] != 0)) encapsulation_load [i] = encapsulation_load [i] - offset_count;
        }
        
        if (std_my_VR=="SQ"){
          encapsulation_level = encapsulation_level + 1;
          encapsulation_load[encapsulation_level] = my_VL;
        }
      }
      
    }
    event++;
  }
  
  
  if (full_info)
    return(Rcpp::List::create(
        Rcpp::Named("error") = "",
        Rcpp::Named("VRness") = std::to_string (decoded_VRness),
        Rcpp::Named("db") = Rcpp::DataFrame::create(
          Rcpp::Named("tag") = vTAG, Rcpp::Named("VR") = vVR, 
          Rcpp::Named("endian") = vendian, Rcpp::Named("start") = vstart,
          Rcpp::Named("stop") = vstop,
          Rcpp::Named("encaps.load") =  vencapsload,
          Rcpp::Named("load.start") =  vloadstart,
          Rcpp::Named("load.stop") =  vloadstop,
          Rcpp::Named("tag.start") =  vtagstart)));
  
  return(Rcpp::List::create(
      Rcpp::Named("error") = "",
      Rcpp::Named("VRness") = std::to_string (decoded_VRness),
      Rcpp::Named("db") = Rcpp::DataFrame::create(
        Rcpp::Named("tag") = vTAG, Rcpp::Named("VR") = vVR, 
        Rcpp::Named("endian") = vendian, Rcpp::Named("start") = vstart,
        Rcpp::Named("stop") = vstop)));
}
