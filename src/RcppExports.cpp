// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// mult_
arma::mat mult_(const arma::mat& X, const arma::mat& Y);
RcppExport SEXP _HOFS_mult_(SEXP XSEXP, SEXP YSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Y(YSEXP);
    rcpp_result_gen = Rcpp::wrap(mult_(X, Y));
    return rcpp_result_gen;
END_RCPP
}
// array_mult
List array_mult(const arma::mat& X, const arma::mat& Y);
RcppExport SEXP _HOFS_array_mult(SEXP XSEXP, SEXP YSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Y(YSEXP);
    rcpp_result_gen = Rcpp::wrap(array_mult(X, Y));
    return rcpp_result_gen;
END_RCPP
}
// Reduce_cpp
NumericMatrix Reduce_cpp(List x);
RcppExport SEXP _HOFS_Reduce_cpp(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(Reduce_cpp(x));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_HOFS_mult_", (DL_FUNC) &_HOFS_mult_, 2},
    {"_HOFS_array_mult", (DL_FUNC) &_HOFS_array_mult, 2},
    {"_HOFS_Reduce_cpp", (DL_FUNC) &_HOFS_Reduce_cpp, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_HOFS(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
