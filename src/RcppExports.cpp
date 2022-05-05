// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// EpiInvertC
List EpiInvertC(NumericVector i_original0, String last_incidence_date, CharacterVector festive_days, NumericVector si_distr0, int shift_si_distr, int max_time_interval, double mean_si, double sd_si, double shift_si, double Rt_regularization_weight, double seasonality_regularization_weight);
RcppExport SEXP _EpiInvertNew_EpiInvertC(SEXP i_original0SEXP, SEXP last_incidence_dateSEXP, SEXP festive_daysSEXP, SEXP si_distr0SEXP, SEXP shift_si_distrSEXP, SEXP max_time_intervalSEXP, SEXP mean_siSEXP, SEXP sd_siSEXP, SEXP shift_siSEXP, SEXP Rt_regularization_weightSEXP, SEXP seasonality_regularization_weightSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type i_original0(i_original0SEXP);
    Rcpp::traits::input_parameter< String >::type last_incidence_date(last_incidence_dateSEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type festive_days(festive_daysSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type si_distr0(si_distr0SEXP);
    Rcpp::traits::input_parameter< int >::type shift_si_distr(shift_si_distrSEXP);
    Rcpp::traits::input_parameter< int >::type max_time_interval(max_time_intervalSEXP);
    Rcpp::traits::input_parameter< double >::type mean_si(mean_siSEXP);
    Rcpp::traits::input_parameter< double >::type sd_si(sd_siSEXP);
    Rcpp::traits::input_parameter< double >::type shift_si(shift_siSEXP);
    Rcpp::traits::input_parameter< double >::type Rt_regularization_weight(Rt_regularization_weightSEXP);
    Rcpp::traits::input_parameter< double >::type seasonality_regularization_weight(seasonality_regularization_weightSEXP);
    rcpp_result_gen = Rcpp::wrap(EpiInvertC(i_original0, last_incidence_date, festive_days, si_distr0, shift_si_distr, max_time_interval, mean_si, sd_si, shift_si, Rt_regularization_weight, seasonality_regularization_weight));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_EpiInvertNew_EpiInvertC", (DL_FUNC) &_EpiInvertNew_EpiInvertC, 11},
    {NULL, NULL, 0}
};

RcppExport void R_init_EpiInvertNew(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
