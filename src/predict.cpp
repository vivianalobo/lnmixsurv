// -*- mode: C++; c-indent-level: 2; c-basic-offset: 2; indent-tabs-mode: nil; -*-

#include <RcppArmadillo.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>

using namespace Rcpp;

// Functions used to predict EM survival
double sob_lognormal_em(const double& t, const double& m, const double& sigma) {
  return R::pnorm((m -log(t))/sigma, 0, 1, true, false);
}

double sob_lognormal_em_mix(const double& t, const arma::rowvec& m, const arma::vec& sigma, const arma::vec& eta) {
  double res = 0.0;
  for (int i = 0; i < m.n_elem; i++) {
    res += eta(i) * sob_lognormal_em(t, m(i), sigma(i));
  }
  return res;
}

double hazard_lognormal_em_mix(const double& t, const arma::rowvec& m, const arma::vec& sigma, const arma::vec& eta) {
  double sob_mix = sob_lognormal_em_mix(t, m, sigma, eta);
  double dlnorm_mix = 0;

  for (int i = 0; i < m.n_elem; i++) {
    dlnorm_mix += eta(i) * R::dlnorm(t, m(i), sigma(i), false);
  }
  return dlnorm_mix/sob_mix;
}

// [[Rcpp::export]]
arma::vec predict_survival_em_cpp(const arma::vec& t, const arma::mat& m, const arma::vec& sigma, const arma::vec& eta, const int& r) {
  int n = t.n_elem;
  arma::vec out(t.n_elem);

  for(int i = 0; i < n; i++) {
    out(i) = sob_lognormal_em_mix(t(i), m.row(r - 1), sigma, eta);
  }
  return out;
}

// [[Rcpp::export]]
arma::vec predict_hazard_em_cpp(const arma::vec& t, const arma::mat& m, const arma::vec& sigma, const arma::vec& eta, const int& r) {
  int n = t.n_elem;
  arma::vec out(t.n_elem);

  for(int i = 0; i < n; i++) {
    out(i) = hazard_lognormal_em_mix(t(i), m.row(r - 1), sigma, eta);
  }
  return out;
}