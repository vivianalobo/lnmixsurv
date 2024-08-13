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
arma::cube lognormal_mixture_gibbs(const int& Niter, const int& em_iter, const int& G, const arma::vec& t, const arma::ivec& delta, const arma::mat& X, const arma::vec& starting_seed, const bool& show_output, const int& n_chains, const bool& use_W, const bool& better_initial_values, const int& N_em, const int& Niter_em, const bool& fast_groups);
RcppExport SEXP _lnmixsurv_lognormal_mixture_gibbs(SEXP NiterSEXP, SEXP em_iterSEXP, SEXP GSEXP, SEXP tSEXP, SEXP deltaSEXP, SEXP XSEXP, SEXP starting_seedSEXP, SEXP show_outputSEXP, SEXP n_chainsSEXP, SEXP use_WSEXP, SEXP better_initial_valuesSEXP, SEXP N_emSEXP, SEXP Niter_emSEXP, SEXP fast_groupsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int& >::type Niter(NiterSEXP);
    Rcpp::traits::input_parameter< const int& >::type em_iter(em_iterSEXP);
    Rcpp::traits::input_parameter< const int& >::type G(GSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type t(tSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type delta(deltaSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type starting_seed(starting_seedSEXP);
    Rcpp::traits::input_parameter< const bool& >::type show_output(show_outputSEXP);
    Rcpp::traits::input_parameter< const int& >::type n_chains(n_chainsSEXP);
    Rcpp::traits::input_parameter< const bool& >::type use_W(use_WSEXP);
    Rcpp::traits::input_parameter< const bool& >::type better_initial_values(better_initial_valuesSEXP);
    Rcpp::traits::input_parameter< const int& >::type N_em(N_emSEXP);
    Rcpp::traits::input_parameter< const int& >::type Niter_em(Niter_emSEXP);
    Rcpp::traits::input_parameter< const bool& >::type fast_groups(fast_groupsSEXP);
    rcpp_result_gen = Rcpp::wrap(lognormal_mixture_gibbs(Niter, em_iter, G, t, delta, X, starting_seed, show_output, n_chains, use_W, better_initial_values, N_em, Niter_em, fast_groups));
    return rcpp_result_gen;
END_RCPP
}
// lognormal_mixture_em_implementation
arma::mat lognormal_mixture_em_implementation(const int& Niter, const int& G, const arma::vec& t, const arma::ivec& delta, const arma::mat& X, long long int starting_seed, const bool& better_initial_values, const int& N_em, const int& Niter_em, const bool& show_output);
RcppExport SEXP _lnmixsurv_lognormal_mixture_em_implementation(SEXP NiterSEXP, SEXP GSEXP, SEXP tSEXP, SEXP deltaSEXP, SEXP XSEXP, SEXP starting_seedSEXP, SEXP better_initial_valuesSEXP, SEXP N_emSEXP, SEXP Niter_emSEXP, SEXP show_outputSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int& >::type Niter(NiterSEXP);
    Rcpp::traits::input_parameter< const int& >::type G(GSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type t(tSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type delta(deltaSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< long long int >::type starting_seed(starting_seedSEXP);
    Rcpp::traits::input_parameter< const bool& >::type better_initial_values(better_initial_valuesSEXP);
    Rcpp::traits::input_parameter< const int& >::type N_em(N_emSEXP);
    Rcpp::traits::input_parameter< const int& >::type Niter_em(Niter_emSEXP);
    Rcpp::traits::input_parameter< const bool& >::type show_output(show_outputSEXP);
    rcpp_result_gen = Rcpp::wrap(lognormal_mixture_em_implementation(Niter, G, t, delta, X, starting_seed, better_initial_values, N_em, Niter_em, show_output));
    return rcpp_result_gen;
END_RCPP
}
// predict_survival_em_cpp
arma::vec predict_survival_em_cpp(const arma::vec& t, const arma::mat& m, const arma::vec& sigma, const arma::vec& eta, const int& r);
RcppExport SEXP _lnmixsurv_predict_survival_em_cpp(SEXP tSEXP, SEXP mSEXP, SEXP sigmaSEXP, SEXP etaSEXP, SEXP rSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type t(tSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type m(mSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< const int& >::type r(rSEXP);
    rcpp_result_gen = Rcpp::wrap(predict_survival_em_cpp(t, m, sigma, eta, r));
    return rcpp_result_gen;
END_RCPP
}
// predict_hazard_em_cpp
arma::vec predict_hazard_em_cpp(const arma::vec& t, const arma::mat& m, const arma::vec& sigma, const arma::vec& eta, const int& r);
RcppExport SEXP _lnmixsurv_predict_hazard_em_cpp(SEXP tSEXP, SEXP mSEXP, SEXP sigmaSEXP, SEXP etaSEXP, SEXP rSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type t(tSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type m(mSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< const int& >::type r(rSEXP);
    rcpp_result_gen = Rcpp::wrap(predict_hazard_em_cpp(t, m, sigma, eta, r));
    return rcpp_result_gen;
END_RCPP
}
// predict_survival_gibbs_cpp
arma::mat predict_survival_gibbs_cpp(const arma::vec& eval_time, const arma::rowvec& predictors, const arma::field<arma::mat>& beta_start, const arma::mat sigma_start, const arma::mat eta_start, const bool& interval, const double& level);
RcppExport SEXP _lnmixsurv_predict_survival_gibbs_cpp(SEXP eval_timeSEXP, SEXP predictorsSEXP, SEXP beta_startSEXP, SEXP sigma_startSEXP, SEXP eta_startSEXP, SEXP intervalSEXP, SEXP levelSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type eval_time(eval_timeSEXP);
    Rcpp::traits::input_parameter< const arma::rowvec& >::type predictors(predictorsSEXP);
    Rcpp::traits::input_parameter< const arma::field<arma::mat>& >::type beta_start(beta_startSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type sigma_start(sigma_startSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type eta_start(eta_startSEXP);
    Rcpp::traits::input_parameter< const bool& >::type interval(intervalSEXP);
    Rcpp::traits::input_parameter< const double& >::type level(levelSEXP);
    rcpp_result_gen = Rcpp::wrap(predict_survival_gibbs_cpp(eval_time, predictors, beta_start, sigma_start, eta_start, interval, level));
    return rcpp_result_gen;
END_RCPP
}
// predict_hazard_gibbs_cpp
arma::mat predict_hazard_gibbs_cpp(const arma::vec& eval_time, const arma::rowvec& predictors, const arma::field<arma::mat>& beta_start, const arma::mat sigma_start, const arma::mat eta_start, const bool& interval, const double& level);
RcppExport SEXP _lnmixsurv_predict_hazard_gibbs_cpp(SEXP eval_timeSEXP, SEXP predictorsSEXP, SEXP beta_startSEXP, SEXP sigma_startSEXP, SEXP eta_startSEXP, SEXP intervalSEXP, SEXP levelSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type eval_time(eval_timeSEXP);
    Rcpp::traits::input_parameter< const arma::rowvec& >::type predictors(predictorsSEXP);
    Rcpp::traits::input_parameter< const arma::field<arma::mat>& >::type beta_start(beta_startSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type sigma_start(sigma_startSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type eta_start(eta_startSEXP);
    Rcpp::traits::input_parameter< const bool& >::type interval(intervalSEXP);
    Rcpp::traits::input_parameter< const double& >::type level(levelSEXP);
    rcpp_result_gen = Rcpp::wrap(predict_hazard_gibbs_cpp(eval_time, predictors, beta_start, sigma_start, eta_start, interval, level));
    return rcpp_result_gen;
END_RCPP
}
// simulate_y
arma::vec simulate_y(const arma::mat& X, const arma::mat& beta, const arma::vec& phi, const arma::ivec& delta, const arma::ivec& groups, long long int starting_seed);
RcppExport SEXP _lnmixsurv_simulate_y(SEXP XSEXP, SEXP betaSEXP, SEXP phiSEXP, SEXP deltaSEXP, SEXP groupsSEXP, SEXP starting_seedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type phi(phiSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type delta(deltaSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type groups(groupsSEXP);
    Rcpp::traits::input_parameter< long long int >::type starting_seed(starting_seedSEXP);
    rcpp_result_gen = Rcpp::wrap(simulate_y(X, beta, phi, delta, groups, starting_seed));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_lnmixsurv_lognormal_mixture_gibbs", (DL_FUNC) &_lnmixsurv_lognormal_mixture_gibbs, 14},
    {"_lnmixsurv_lognormal_mixture_em_implementation", (DL_FUNC) &_lnmixsurv_lognormal_mixture_em_implementation, 10},
    {"_lnmixsurv_predict_survival_em_cpp", (DL_FUNC) &_lnmixsurv_predict_survival_em_cpp, 5},
    {"_lnmixsurv_predict_hazard_em_cpp", (DL_FUNC) &_lnmixsurv_predict_hazard_em_cpp, 5},
    {"_lnmixsurv_predict_survival_gibbs_cpp", (DL_FUNC) &_lnmixsurv_predict_survival_gibbs_cpp, 7},
    {"_lnmixsurv_predict_hazard_gibbs_cpp", (DL_FUNC) &_lnmixsurv_predict_hazard_gibbs_cpp, 7},
    {"_lnmixsurv_simulate_y", (DL_FUNC) &_lnmixsurv_simulate_y, 6},
    {NULL, NULL, 0}
};

RcppExport void R_init_lnmixsurv(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
