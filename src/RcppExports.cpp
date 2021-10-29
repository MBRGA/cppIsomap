// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// cppIsomap
List cppIsomap(NumericMatrix data, IntegerVector dims, Nullable<IntegerVector> k, bool mod, bool verbose);
RcppExport SEXP _cppIsomap_cppIsomap(SEXP dataSEXP, SEXP dimsSEXP, SEXP kSEXP, SEXP modSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type data(dataSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type dims(dimsSEXP);
    Rcpp::traits::input_parameter< Nullable<IntegerVector> >::type k(kSEXP);
    Rcpp::traits::input_parameter< bool >::type mod(modSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(cppIsomap(data, dims, k, mod, verbose));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_cppIsomap_cppIsomap", (DL_FUNC) &_cppIsomap_cppIsomap, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_cppIsomap(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
