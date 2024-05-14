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

// lognormal_mixture_gibbs
arma::cube lognormal_mixture_gibbs(const int& Niter, const int& em_iter, const int& G, const arma::vec& exp_y, const arma::ivec& delta, const arma::mat& X, const double& a, arma::Col<long long int> starting_seed, const bool& show_output, const int& n_cores, const int& n_chains, const bool& force_num_cores, const bool& sparse, const bool& use_W, const bool& better_initial_values, const int& N_em, const int& Niter_em);
RcppExport SEXP _lnmixsurv_lognormal_mixture_gibbs(SEXP NiterSEXP, SEXP em_iterSEXP, SEXP GSEXP, SEXP exp_ySEXP, SEXP deltaSEXP, SEXP XSEXP, SEXP aSEXP, SEXP starting_seedSEXP, SEXP show_outputSEXP, SEXP n_coresSEXP, SEXP n_chainsSEXP, SEXP force_num_coresSEXP, SEXP sparseSEXP, SEXP use_WSEXP, SEXP better_initial_valuesSEXP, SEXP N_emSEXP, SEXP Niter_emSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int& >::type Niter(NiterSEXP);
    Rcpp::traits::input_parameter< const int& >::type em_iter(em_iterSEXP);
    Rcpp::traits::input_parameter< const int& >::type G(GSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type exp_y(exp_ySEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type delta(deltaSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const double& >::type a(aSEXP);
    Rcpp::traits::input_parameter< arma::Col<long long int> >::type starting_seed(starting_seedSEXP);
    Rcpp::traits::input_parameter< const bool& >::type show_output(show_outputSEXP);
    Rcpp::traits::input_parameter< const int& >::type n_cores(n_coresSEXP);
    Rcpp::traits::input_parameter< const int& >::type n_chains(n_chainsSEXP);
    Rcpp::traits::input_parameter< const bool& >::type force_num_cores(force_num_coresSEXP);
    Rcpp::traits::input_parameter< const bool& >::type sparse(sparseSEXP);
    Rcpp::traits::input_parameter< const bool& >::type use_W(use_WSEXP);
    Rcpp::traits::input_parameter< const bool& >::type better_initial_values(better_initial_valuesSEXP);
    Rcpp::traits::input_parameter< const int& >::type N_em(N_emSEXP);
    Rcpp::traits::input_parameter< const int& >::type Niter_em(Niter_emSEXP);
    rcpp_result_gen = Rcpp::wrap(lognormal_mixture_gibbs(Niter, em_iter, G, exp_y, delta, X, a, starting_seed, show_output, n_cores, n_chains, force_num_cores, sparse, use_W, better_initial_values, N_em, Niter_em));
    return rcpp_result_gen;
END_RCPP
}
// lognormal_mixture_em_implementation
arma::mat lognormal_mixture_em_implementation(const int& Niter, const int& G, const arma::vec& t, const arma::ivec& delta, const arma::mat& X, long long int starting_seed, const bool& sparse, const bool& better_initial_values, const int& N_em, const int& Niter_em);
RcppExport SEXP _lnmixsurv_lognormal_mixture_em_implementation(SEXP NiterSEXP, SEXP GSEXP, SEXP tSEXP, SEXP deltaSEXP, SEXP XSEXP, SEXP starting_seedSEXP, SEXP sparseSEXP, SEXP better_initial_valuesSEXP, SEXP N_emSEXP, SEXP Niter_emSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int& >::type Niter(NiterSEXP);
    Rcpp::traits::input_parameter< const int& >::type G(GSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type t(tSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type delta(deltaSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< long long int >::type starting_seed(starting_seedSEXP);
    Rcpp::traits::input_parameter< const bool& >::type sparse(sparseSEXP);
    Rcpp::traits::input_parameter< const bool& >::type better_initial_values(better_initial_valuesSEXP);
    Rcpp::traits::input_parameter< const int& >::type N_em(N_emSEXP);
    Rcpp::traits::input_parameter< const int& >::type Niter_em(Niter_emSEXP);
    rcpp_result_gen = Rcpp::wrap(lognormal_mixture_em_implementation(Niter, G, t, delta, X, starting_seed, sparse, better_initial_values, N_em, Niter_em));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_lnmixsurv_lognormal_mixture_gibbs", (DL_FUNC) &_lnmixsurv_lognormal_mixture_gibbs, 17},
    {"_lnmixsurv_lognormal_mixture_em_implementation", (DL_FUNC) &_lnmixsurv_lognormal_mixture_em_implementation, 10},
    {NULL, NULL, 0}
};

RcppExport void R_init_lnmixsurv(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
