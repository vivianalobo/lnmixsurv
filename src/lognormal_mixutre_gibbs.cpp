// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// [[Rcpp::plugins("cpp11")]]

#include <RcppArmadillo.h>
#include <mvnorm.h>
#include "RcppTN.h"
#include <progress.hpp>
#include <progress_bar.hpp>

arma::colvec dnorm_vec(arma::colvec y, arma::colvec mean, double sd, bool log)
{
    int n = y.size();
    arma::colvec ret(n);
    for (int i = 0; i < n; i++)
    {
        ret(i) = R::dnorm(y(i), mean(i), sd, log);
    }
    return ret;
}

// Realiza n amostras de uma bernoulli, cada uma com uma prob
// diferente. prob deve ter tamanho n.
arma::colvec rbernoulli(int n, arma::colvec prob)
{
    arma::colvec ret(n);
    for (int i = 0; i < n; i++)
    {
        ret(i) = R::rbinom(1, prob(i));
    }
    return ret;
}

// rtn1 mas com alguns argumentos especificos vetorizados.
arma::colvec rtn1_vec(int n, arma::colvec mean, double sd, arma::colvec low, double high)
{
    arma::colvec ret(n);
    for (int i = 0; i < n; i++)
    {
        ret(i) = RcppTN::rtn1(mean(i), sd, low(i), high);
    }
    return ret;
}

arma::colvec calcular_prob(arma::colvec log_y, double theta, arma::mat x, arma::colvec t_beta_a, arma::colvec t_beta_b, double sd_a, double sd_b)
{
    arma::colvec mean_a = x * t_beta_a;
    arma::colvec mean_b = x * t_beta_b;
    arma::colvec pr1 = theta * dnorm_vec(log_y, mean_a, sd_a, false);
    arma::colvec pr2 = (1 - theta) * dnorm_vec(log_y, mean_b, sd_b, false);
    return pr1 / (pr1 + pr2);
}

//' Realiza passo da augmentation.
//'
//' O valor de log_y eh alterado no procedimento.
//'
//' @param indicadora_grupo [nx1] vetor de zeros e uns
//' @param grupo [1x1] grupo ao qual sera realizada o augmentation (0 ou 1)
//' @param log_y [nx1]
//' @param c [nx1] uma copia do valor original de log_y.
//' @param x [nxp] matriz de covariaveis
//' @param cens [nx1] vetor indicador de censura (1 = censura, 0 = evento).
void realizar_augmentation(arma::uvec idx, arma::colvec& log_y, arma::colvec c, arma::mat x, arma::uvec cens, arma::colvec t_beta, double sd, arma::colvec indicadora_grupo, int grupo)
{
    arma::uvec g_z = arma::find(indicadora_grupo == grupo && cens == 1);
    int tamanho = g_z.size();

    // augmentation
    arma::colvec g_z_mean = x.rows(g_z) * t_beta;

    log_y.elem(g_z) = rtn1_vec(tamanho, g_z_mean, sd, c.elem(g_z), R_PosInf);
    return;
}

//' Gibbs sampling for a two componentes mixture model
//'
//' @param y [nx1]
//' @param x [nxp]
//' @param delta [nx1]
//'
// [[Rcpp::export]]
Rcpp::List lognormal_mixture_gibbs_cpp(arma::colvec y, arma::mat x, arma::colvec delta, int numero_iteracoes, double valor_inicial_beta = 0)
{

    int numero_observacoes = y.size();
    int numero_covariaveis = x.n_cols;

    arma::colvec log_y = log(y);
    arma::uvec cens = (delta == 0);
    arma::colvec theta(numero_iteracoes);
    arma::colvec phi_a(numero_iteracoes);
    arma::colvec phi_b(numero_iteracoes);
    arma::mat beta_a(numero_iteracoes, numero_covariaveis);
    arma::mat beta_b(numero_iteracoes, numero_covariaveis);
    arma::colvec c = arma::colvec(log_y);

    // valores iniciais
    theta(0) = R::runif(0, 1);
    phi_a(0) = R::runif(0, 1);
    phi_b(0) = R::runif(0, 1);

    beta_a.row(0) = arma::rowvec(numero_covariaveis, arma::fill::value(valor_inicial_beta));
    beta_b.row(0) = arma::rowvec(numero_covariaveis, arma::fill::value(valor_inicial_beta));

    // variaveis auxiliares usadas no for
    arma::colvec mean_a(numero_observacoes);
    arma::colvec mean_b(numero_observacoes);
    double sd_a;
    double sd_b;
    arma::colvec prob(numero_observacoes);
    arma::uvec idxA;
    arma::uvec idxB;
    arma::mat XA;
    arma::mat XB;
    arma::colvec yA;
    arma::colvec yB;
    arma::colvec I(numero_observacoes);
    int A;
    int B;
    arma::colvec c1A;
    double c2A;
    arma::colvec c1B;
    double c2B;

    arma::mat tXA;
    arma::mat tXB;

    arma::mat vA;
    arma::mat vB;

    arma::colvec mA;
    arma::colvec mB;

    arma::colvec aux_beta_a(numero_covariaveis);
    arma::colvec aux_beta_b(numero_covariaveis);

    Progress pro(numero_iteracoes, true);
    for (int it = 1; it <= numero_iteracoes - 1; it++)
    {
        if (Progress::check_abort())
            return Rcpp::List::create(
                Rcpp::_["beta_a"] = beta_a,
                Rcpp::_["beta_b"] = beta_b,
                Rcpp::_["phi_a"] = phi_a,
                Rcpp::_["phi_b"] = phi_b,
                Rcpp::_["theta"] = theta);
        pro.increment();

        aux_beta_a = arma::trans(beta_a.row(it - 1));
        aux_beta_b = arma::trans(beta_b.row(it - 1));

        sd_a = 1 / sqrt(phi_a(it - 1));
        sd_b = 1 / sqrt(phi_b(it - 1));

        // probability
        prob = calcular_prob(log_y, theta(it - 1), x, aux_beta_a, aux_beta_b, sd_a, sd_b);

        // mixture
        I = rbernoulli(numero_observacoes, prob);
        idxA = arma::find(I == 1);
        idxB = arma::find(I == 0);
        XA = x.rows(idxA);
        XB = x.rows(idxB);
        yA = log_y.elem(idxA);
        yB = log_y.elem(idxB);

        realizar_augmentation(idxA, log_y, c, x, cens, aux_beta_a, sd_a, I, 1);
        realizar_augmentation(idxB, log_y, c, x, cens, aux_beta_b, sd_b, I, 0);

        // theta
        A = idxA.size();
        B = idxB.size();
        theta(it) = R::rbeta(1 + A, 1 + B);

        // phi
        c1A = yA - XA * aux_beta_a;
        c2A = arma::dot(c1A, c1A);
        c1B = yB - XB * aux_beta_b;
        c2B = arma::dot(c1B, c1B);

        phi_a(it) = R::rgamma(0.01 + A / 2, 1 / (0.01 + c2A / 2));
        phi_b(it) = R::rgamma(0.01 + B / 2, 1 / (0.01 + c2B / 2));

        tXA = arma::trans(XA);
        vA = arma::inv(tXA * XA, arma::inv_opts::allow_approx);
        mA = vA * tXA * yA;
        tXB = arma::trans(XB);
        vB = arma::inv(tXB * XB, arma::inv_opts::allow_approx);
        mB = vB * tXB * yB;

        beta_a.row(it) = rmvnorm(1, mA, vA / phi_a(it));
        beta_b.row(it) = rmvnorm(1, mB, vB / phi_b(it));
    }
    return Rcpp::List::create(
        Rcpp::_["betaA"] = beta_a,
        Rcpp::_["betaB"] = beta_b,
        Rcpp::_["phiA"] = phi_a,
        Rcpp::_["phiB"] = phi_b,
        Rcpp::_["theta"] = theta);
}
