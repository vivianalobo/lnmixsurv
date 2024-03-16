// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <RcppGSL.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// lognormal_mixture_em
arma::mat lognormal_mixture_em(int Niter, int G, arma::vec y, arma::vec delta, arma::mat X, long long int starting_seed);
RcppExport SEXP _lnmixsurv_lognormal_mixture_em(SEXP NiterSEXP, SEXP GSEXP, SEXP ySEXP, SEXP deltaSEXP, SEXP XSEXP, SEXP starting_seedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type Niter(NiterSEXP);
    Rcpp::traits::input_parameter< int >::type G(GSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::vec >::type delta(deltaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< long long int >::type starting_seed(starting_seedSEXP);
    rcpp_result_gen = Rcpp::wrap(lognormal_mixture_em(Niter, G, y, delta, X, starting_seed));
    return rcpp_result_gen;
END_RCPP
}
// lognormal_mixture_gibbs
arma::cube lognormal_mixture_gibbs(int Niter, int em_iter, int G, arma::vec exp_y, arma::ivec delta, arma::mat X, double a, arma::Col<long long int> starting_seed, bool show_output, int n_cores, int n_chains, bool force_num_cores);
RcppExport SEXP _lnmixsurv_lognormal_mixture_gibbs(SEXP NiterSEXP, SEXP em_iterSEXP, SEXP GSEXP, SEXP exp_ySEXP, SEXP deltaSEXP, SEXP XSEXP, SEXP aSEXP, SEXP starting_seedSEXP, SEXP show_outputSEXP, SEXP n_coresSEXP, SEXP n_chainsSEXP, SEXP force_num_coresSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type Niter(NiterSEXP);
    Rcpp::traits::input_parameter< int >::type em_iter(em_iterSEXP);
    Rcpp::traits::input_parameter< int >::type G(GSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type exp_y(exp_ySEXP);
    Rcpp::traits::input_parameter< arma::ivec >::type delta(deltaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< arma::Col<long long int> >::type starting_seed(starting_seedSEXP);
    Rcpp::traits::input_parameter< bool >::type show_output(show_outputSEXP);
    Rcpp::traits::input_parameter< int >::type n_cores(n_coresSEXP);
    Rcpp::traits::input_parameter< int >::type n_chains(n_chainsSEXP);
    Rcpp::traits::input_parameter< bool >::type force_num_cores(force_num_coresSEXP);
    rcpp_result_gen = Rcpp::wrap(lognormal_mixture_gibbs(Niter, em_iter, G, exp_y, delta, X, a, starting_seed, show_output, n_cores, n_chains, force_num_cores));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_lnmixsurv_lognormal_mixture_em", (DL_FUNC) &_lnmixsurv_lognormal_mixture_em, 6},
    {"_lnmixsurv_lognormal_mixture_gibbs", (DL_FUNC) &_lnmixsurv_lognormal_mixture_gibbs, 12},
    {NULL, NULL, 0}
};

RcppExport void R_init_lnmixsurv(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
