// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// dicombrowser
std::list <std::string> dicombrowser(std::vector <unsigned char> dicomrawdata, DataFrame tagdico, unsigned int nbTAG, std::string stop_tag, unsigned int stop_level, bool full_info, bool verbose);
RcppExport SEXP _espadon_dicombrowser(SEXP dicomrawdataSEXP, SEXP tagdicoSEXP, SEXP nbTAGSEXP, SEXP stop_tagSEXP, SEXP stop_levelSEXP, SEXP full_infoSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector <unsigned char> >::type dicomrawdata(dicomrawdataSEXP);
    Rcpp::traits::input_parameter< DataFrame >::type tagdico(tagdicoSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type nbTAG(nbTAGSEXP);
    Rcpp::traits::input_parameter< std::string >::type stop_tag(stop_tagSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type stop_level(stop_levelSEXP);
    Rcpp::traits::input_parameter< bool >::type full_info(full_infoSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(dicombrowser(dicomrawdata, tagdico, nbTAG, stop_tag, stop_level, full_info, verbose));
    return rcpp_result_gen;
END_RCPP
}
// fansphereC
std::vector <double> fansphereC(double angle);
RcppExport SEXP _espadon_fansphereC(SEXP angleSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type angle(angleSEXP);
    rcpp_result_gen = Rcpp::wrap(fansphereC(angle));
    return rcpp_result_gen;
END_RCPP
}
// fantovoxelC
std::vector <double> fantovoxelC(std::vector <double> p, std::vector  <int> n_ijk, std::vector  <int> k_idx, std::vector  <int> k_loc, std::vector <double> O_ijk, std::vector <double> vol_data, bool att, bool vol_value_flag, double vol_value);
RcppExport SEXP _espadon_fantovoxelC(SEXP pSEXP, SEXP n_ijkSEXP, SEXP k_idxSEXP, SEXP k_locSEXP, SEXP O_ijkSEXP, SEXP vol_dataSEXP, SEXP attSEXP, SEXP vol_value_flagSEXP, SEXP vol_valueSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector <double> >::type p(pSEXP);
    Rcpp::traits::input_parameter< std::vector  <int> >::type n_ijk(n_ijkSEXP);
    Rcpp::traits::input_parameter< std::vector  <int> >::type k_idx(k_idxSEXP);
    Rcpp::traits::input_parameter< std::vector  <int> >::type k_loc(k_locSEXP);
    Rcpp::traits::input_parameter< std::vector <double> >::type O_ijk(O_ijkSEXP);
    Rcpp::traits::input_parameter< std::vector <double> >::type vol_data(vol_dataSEXP);
    Rcpp::traits::input_parameter< bool >::type att(attSEXP);
    Rcpp::traits::input_parameter< bool >::type vol_value_flag(vol_value_flagSEXP);
    Rcpp::traits::input_parameter< double >::type vol_value(vol_valueSEXP);
    rcpp_result_gen = Rcpp::wrap(fantovoxelC(p, n_ijk, k_idx, k_loc, O_ijk, vol_data, att, vol_value_flag, vol_value));
    return rcpp_result_gen;
END_RCPP
}
// gammaindex
std::vector <double> gammaindex(std::vector <double> vol3D, std::vector <double> vol3D_ref, std::vector <int> inspect_idx, std::vector <int> n_ijk, std::vector <double> rel_dxyz, std::vector  <int> ball_i, std::vector  <int> ball_j, std::vector  <int> ball_k, int around_idx, std::vector  <double> distance, double D_norm, bool local, double local_th_pc, double ref_pc);
RcppExport SEXP _espadon_gammaindex(SEXP vol3DSEXP, SEXP vol3D_refSEXP, SEXP inspect_idxSEXP, SEXP n_ijkSEXP, SEXP rel_dxyzSEXP, SEXP ball_iSEXP, SEXP ball_jSEXP, SEXP ball_kSEXP, SEXP around_idxSEXP, SEXP distanceSEXP, SEXP D_normSEXP, SEXP localSEXP, SEXP local_th_pcSEXP, SEXP ref_pcSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector <double> >::type vol3D(vol3DSEXP);
    Rcpp::traits::input_parameter< std::vector <double> >::type vol3D_ref(vol3D_refSEXP);
    Rcpp::traits::input_parameter< std::vector <int> >::type inspect_idx(inspect_idxSEXP);
    Rcpp::traits::input_parameter< std::vector <int> >::type n_ijk(n_ijkSEXP);
    Rcpp::traits::input_parameter< std::vector <double> >::type rel_dxyz(rel_dxyzSEXP);
    Rcpp::traits::input_parameter< std::vector  <int> >::type ball_i(ball_iSEXP);
    Rcpp::traits::input_parameter< std::vector  <int> >::type ball_j(ball_jSEXP);
    Rcpp::traits::input_parameter< std::vector  <int> >::type ball_k(ball_kSEXP);
    Rcpp::traits::input_parameter< int >::type around_idx(around_idxSEXP);
    Rcpp::traits::input_parameter< std::vector  <double> >::type distance(distanceSEXP);
    Rcpp::traits::input_parameter< double >::type D_norm(D_normSEXP);
    Rcpp::traits::input_parameter< bool >::type local(localSEXP);
    Rcpp::traits::input_parameter< double >::type local_th_pc(local_th_pcSEXP);
    Rcpp::traits::input_parameter< double >::type ref_pc(ref_pcSEXP);
    rcpp_result_gen = Rcpp::wrap(gammaindex(vol3D, vol3D_ref, inspect_idx, n_ijk, rel_dxyz, ball_i, ball_j, ball_k, around_idx, distance, D_norm, local, local_th_pc, ref_pc));
    return rcpp_result_gen;
END_RCPP
}
// getijktfromindexC
std::vector <int> getijktfromindexC(std::vector <int> index, std::vector  <int> k_idx, std::vector  <int> n_ijk);
RcppExport SEXP _espadon_getijktfromindexC(SEXP indexSEXP, SEXP k_idxSEXP, SEXP n_ijkSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector <int> >::type index(indexSEXP);
    Rcpp::traits::input_parameter< std::vector  <int> >::type k_idx(k_idxSEXP);
    Rcpp::traits::input_parameter< std::vector  <int> >::type n_ijk(n_ijkSEXP);
    rcpp_result_gen = Rcpp::wrap(getijktfromindexC(index, k_idx, n_ijk));
    return rcpp_result_gen;
END_RCPP
}
// getvaluefromijkC
std::vector <double> getvaluefromijkC(std::vector <double> vol3D, bool interpolate, std::vector <double> i, std::vector <double> j, std::vector <double> k, std::vector  <int> k_idx, std::vector  <int> k_loc, std::vector  <int> n_ijk);
RcppExport SEXP _espadon_getvaluefromijkC(SEXP vol3DSEXP, SEXP interpolateSEXP, SEXP iSEXP, SEXP jSEXP, SEXP kSEXP, SEXP k_idxSEXP, SEXP k_locSEXP, SEXP n_ijkSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector <double> >::type vol3D(vol3DSEXP);
    Rcpp::traits::input_parameter< bool >::type interpolate(interpolateSEXP);
    Rcpp::traits::input_parameter< std::vector <double> >::type i(iSEXP);
    Rcpp::traits::input_parameter< std::vector <double> >::type j(jSEXP);
    Rcpp::traits::input_parameter< std::vector <double> >::type k(kSEXP);
    Rcpp::traits::input_parameter< std::vector  <int> >::type k_idx(k_idxSEXP);
    Rcpp::traits::input_parameter< std::vector  <int> >::type k_loc(k_locSEXP);
    Rcpp::traits::input_parameter< std::vector  <int> >::type n_ijk(n_ijkSEXP);
    rcpp_result_gen = Rcpp::wrap(getvaluefromijkC(vol3D, interpolate, i, j, k, k_idx, k_loc, n_ijk));
    return rcpp_result_gen;
END_RCPP
}
// labelbrowser
std::vector <unsigned int> labelbrowser(std::vector <bool> vol3D, std::vector  <unsigned int> n_ijk);
RcppExport SEXP _espadon_labelbrowser(SEXP vol3DSEXP, SEXP n_ijkSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector <bool> >::type vol3D(vol3DSEXP);
    Rcpp::traits::input_parameter< std::vector  <unsigned int> >::type n_ijk(n_ijkSEXP);
    rcpp_result_gen = Rcpp::wrap(labelbrowser(vol3D, n_ijk));
    return rcpp_result_gen;
END_RCPP
}
// mean_voxC
std::vector <double> mean_voxC(std::vector <double> vol3D, std::vector <int> index, std::vector <int> index_list, std::vector <double> value_list, std::vector <double> value_att_list);
RcppExport SEXP _espadon_mean_voxC(SEXP vol3DSEXP, SEXP indexSEXP, SEXP index_listSEXP, SEXP value_listSEXP, SEXP value_att_listSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector <double> >::type vol3D(vol3DSEXP);
    Rcpp::traits::input_parameter< std::vector <int> >::type index(indexSEXP);
    Rcpp::traits::input_parameter< std::vector <int> >::type index_list(index_listSEXP);
    Rcpp::traits::input_parameter< std::vector <double> >::type value_list(value_listSEXP);
    Rcpp::traits::input_parameter< std::vector <double> >::type value_att_list(value_att_listSEXP);
    rcpp_result_gen = Rcpp::wrap(mean_voxC(vol3D, index, index_list, value_list, value_att_list));
    return rcpp_result_gen;
END_RCPP
}
// medianfilterC
std::vector <double> medianfilterC(std::vector <double> vol3D, std::vector <int> n_ijk, std::vector <long> analyse_idx_vect, std::vector  <int> ball_i, std::vector  <int> ball_j, std::vector  <int> ball_k);
RcppExport SEXP _espadon_medianfilterC(SEXP vol3DSEXP, SEXP n_ijkSEXP, SEXP analyse_idx_vectSEXP, SEXP ball_iSEXP, SEXP ball_jSEXP, SEXP ball_kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector <double> >::type vol3D(vol3DSEXP);
    Rcpp::traits::input_parameter< std::vector <int> >::type n_ijk(n_ijkSEXP);
    Rcpp::traits::input_parameter< std::vector <long> >::type analyse_idx_vect(analyse_idx_vectSEXP);
    Rcpp::traits::input_parameter< std::vector  <int> >::type ball_i(ball_iSEXP);
    Rcpp::traits::input_parameter< std::vector  <int> >::type ball_j(ball_jSEXP);
    Rcpp::traits::input_parameter< std::vector  <int> >::type ball_k(ball_kSEXP);
    rcpp_result_gen = Rcpp::wrap(medianfilterC(vol3D, n_ijk, analyse_idx_vect, ball_i, ball_j, ball_k));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_espadon_dicombrowser", (DL_FUNC) &_espadon_dicombrowser, 7},
    {"_espadon_fansphereC", (DL_FUNC) &_espadon_fansphereC, 1},
    {"_espadon_fantovoxelC", (DL_FUNC) &_espadon_fantovoxelC, 9},
    {"_espadon_gammaindex", (DL_FUNC) &_espadon_gammaindex, 14},
    {"_espadon_getijktfromindexC", (DL_FUNC) &_espadon_getijktfromindexC, 3},
    {"_espadon_getvaluefromijkC", (DL_FUNC) &_espadon_getvaluefromijkC, 8},
    {"_espadon_labelbrowser", (DL_FUNC) &_espadon_labelbrowser, 2},
    {"_espadon_mean_voxC", (DL_FUNC) &_espadon_mean_voxC, 5},
    {"_espadon_medianfilterC", (DL_FUNC) &_espadon_medianfilterC, 6},
    {NULL, NULL, 0}
};

RcppExport void R_init_espadon(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
