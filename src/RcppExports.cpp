// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// wdm_cpp
double wdm_cpp(const std::vector<double>& x, const std::vector<double>& y, std::string method, const std::vector<double>& weights, bool remove_missing);
RcppExport SEXP _wdm_wdm_cpp(SEXP xSEXP, SEXP ySEXP, SEXP methodSEXP, SEXP weightsSEXP, SEXP remove_missingSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<double>& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type y(ySEXP);
    Rcpp::traits::input_parameter< std::string >::type method(methodSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< bool >::type remove_missing(remove_missingSEXP);
    rcpp_result_gen = Rcpp::wrap(wdm_cpp(x, y, method, weights, remove_missing));
    return rcpp_result_gen;
END_RCPP
}
// wdm_mat_cpp
Rcpp::NumericMatrix wdm_mat_cpp(const Rcpp::NumericMatrix& x, std::string method, const std::vector<double>& weights);
RcppExport SEXP _wdm_wdm_mat_cpp(SEXP xSEXP, SEXP methodSEXP, SEXP weightsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type x(xSEXP);
    Rcpp::traits::input_parameter< std::string >::type method(methodSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type weights(weightsSEXP);
    rcpp_result_gen = Rcpp::wrap(wdm_mat_cpp(x, method, weights));
    return rcpp_result_gen;
END_RCPP
}
// indep_test_cpp
Rcpp::List indep_test_cpp(const std::vector<double>& x, const std::vector<double>& y, std::string method, const std::vector<double>& weights, bool remove_missing, std::string alternative);
RcppExport SEXP _wdm_indep_test_cpp(SEXP xSEXP, SEXP ySEXP, SEXP methodSEXP, SEXP weightsSEXP, SEXP remove_missingSEXP, SEXP alternativeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<double>& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type y(ySEXP);
    Rcpp::traits::input_parameter< std::string >::type method(methodSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< bool >::type remove_missing(remove_missingSEXP);
    Rcpp::traits::input_parameter< std::string >::type alternative(alternativeSEXP);
    rcpp_result_gen = Rcpp::wrap(indep_test_cpp(x, y, method, weights, remove_missing, alternative));
    return rcpp_result_gen;
END_RCPP
}
// test
void test();
RcppExport SEXP _wdm_test() {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    test();
    return R_NilValue;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_wdm_wdm_cpp", (DL_FUNC) &_wdm_wdm_cpp, 5},
    {"_wdm_wdm_mat_cpp", (DL_FUNC) &_wdm_wdm_mat_cpp, 3},
    {"_wdm_indep_test_cpp", (DL_FUNC) &_wdm_indep_test_cpp, 6},
    {"_wdm_test", (DL_FUNC) &_wdm_test, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_wdm(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
