// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "../inst/include/action_muv.h"
#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// run_ACTION_muV
List run_ACTION_muV(const List& S, int k_min, int k_max, vec alpha, double lambda, int AA_iters, int Opt_iters, int numThreads);
RcppExport SEXP _ACTIONmv_run_ACTION_muV(SEXP SSEXP, SEXP k_minSEXP, SEXP k_maxSEXP, SEXP alphaSEXP, SEXP lambdaSEXP, SEXP AA_itersSEXP, SEXP Opt_itersSEXP, SEXP numThreadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const List& >::type S(SSEXP);
    Rcpp::traits::input_parameter< int >::type k_min(k_minSEXP);
    Rcpp::traits::input_parameter< int >::type k_max(k_maxSEXP);
    Rcpp::traits::input_parameter< vec >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< int >::type AA_iters(AA_itersSEXP);
    Rcpp::traits::input_parameter< int >::type Opt_iters(Opt_itersSEXP);
    Rcpp::traits::input_parameter< int >::type numThreads(numThreadsSEXP);
    rcpp_result_gen = Rcpp::wrap(run_ACTION_muV(S, k_min, k_max, alpha, lambda, AA_iters, Opt_iters, numThreads));
    return rcpp_result_gen;
END_RCPP
}
// run_test
List run_test(mat H1, mat H2);
RcppExport SEXP _ACTIONmv_run_test(SEXP H1SEXP, SEXP H2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< mat >::type H1(H1SEXP);
    Rcpp::traits::input_parameter< mat >::type H2(H2SEXP);
    rcpp_result_gen = Rcpp::wrap(run_test(H1, H2));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_ACTIONmv_run_ACTION_muV", (DL_FUNC) &_ACTIONmv_run_ACTION_muV, 8},
    {"_ACTIONmv_run_test", (DL_FUNC) &_ACTIONmv_run_test, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_ACTIONmv(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
