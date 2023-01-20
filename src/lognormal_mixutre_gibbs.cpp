// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// [[Rcpp::plugins("cpp11")]]

#include <RcppArmadillo.h>
#include <mvnorm.h>
#include "RcppTN.h"
#include <progress.hpp>
#include <progress_bar.hpp>

#include <RcppClock.h>

arma::colvec dnorm_vec(arma::colvec y, arma::colvec mean, double sd, bool log)
{
    int n = y.size();
    arma::colvec sd_vec(n, arma::fill::value(sd));
    arma::colvec ret = arma::normpdf(y, mean, sd_vec);
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
void realizar_augmentation(arma::uvec idx, arma::colvec &log_y, arma::colvec c, arma::mat x, arma::uvec cens, arma::colvec t_beta, double sd, arma::colvec indicadora_grupo, int grupo)
{
    arma::uvec g_z = arma::find(indicadora_grupo == grupo && cens == 1);
    int tamanho = g_z.size();

    // augmentation
    arma::colvec g_z_mean = x.rows(g_z) * t_beta;

    log_y.elem(g_z) = rtn1_vec(tamanho, g_z_mean, sd, c.elem(g_z), R_PosInf);
    return;
}

double gerar_phi(int n, arma::colvec y, arma::mat x, arma::colvec t_beta)
{
    arma::colvec c1A = y - x * t_beta;
    double c2A = arma::dot(c1A, c1A);
    return R::rgamma(0.01 + n / 2, 1 / (0.01 + c2A / 2));
}

arma::colvec gerar_beta(arma::colvec y, arma::mat x, double phi)
{
    arma::mat tx = arma::trans(x);
    arma::mat v = arma::inv_sympd(tx * x, arma::inv_opts::allow_approx);
    arma::colvec m = v * tx * y;

    return arma::trans(rmvnorm(1, m, v / phi));
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
    arma::mat beta_a(numero_covariaveis, numero_iteracoes);
    arma::mat beta_b(numero_covariaveis, numero_iteracoes);
    arma::colvec c = arma::colvec(log_y);

    // valores iniciais
    theta(0) = R::runif(0, 1);
    phi_a(0) = R::runif(0, 1);
    phi_b(0) = R::runif(0, 1);

    beta_a.col(0) = arma::colvec(numero_covariaveis, arma::fill::value(valor_inicial_beta));
    beta_b.col(0) = arma::colvec(numero_covariaveis, arma::fill::value(valor_inicial_beta));

    // variaveis auxiliares usadas no for
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

    arma::colvec aux_beta_a(numero_covariaveis);
    arma::colvec aux_beta_b(numero_covariaveis);

    // profiling
    // Rcpp::Clock clock;

    Progress pro(numero_iteracoes, true);
    // // clock.tick("algoritmo_completo");
    for (int it = 1; it <= numero_iteracoes - 1; it++)
    {
        // clock.tick("check_abort");
        if (Progress::check_abort())
        {
            // clock.stop("run_times");
            return Rcpp::List::create(
                Rcpp::_["beta_a"] = beta_a,
                Rcpp::_["beta_b"] = beta_b,
                Rcpp::_["phi_a"] = phi_a,
                Rcpp::_["phi_b"] = phi_b,
                Rcpp::_["theta"] = theta);
        }

        // clock.tock("check_abort");
        // clock.tick("increment_pb");
        pro.increment();
        // clock.tock("increment_pb");

        // clock.tick("criar_aux_beta");
        aux_beta_a = beta_a.col(it - 1);
        aux_beta_b = beta_b.col(it - 1);
        // clock.tock("criar_aux_beta");

        // clock.tick("criar_sd");
        sd_a = 1 / sqrt(phi_a(it - 1));
        sd_b = 1 / sqrt(phi_b(it - 1));
        // clock.tock("criar_sd");

        // probability
        // clock.tick("calcular_prob");
        prob = calcular_prob(log_y, theta(it - 1), x, aux_beta_a, aux_beta_b, sd_a, sd_b);
        // clock.tock("calcular_prob");
        // mixture
        // clock.tick("gerar_bernoulli");
        I = rbernoulli(numero_observacoes, prob);
        // clock.tock("gerar_bernoulli");

        // clock.tick("separar_grupos");
        idxA = arma::find(I == 1);
        idxB = arma::find(I == 0);
        // clock.tock("separar_grupos");

        // clock.tick("realizar_augmentation");
        realizar_augmentation(idxA, log_y, c, x, cens, aux_beta_a, sd_a, I, 1);
        realizar_augmentation(idxB, log_y, c, x, cens, aux_beta_b, sd_b, I, 0);
        // clock.tock("realizar_augmentation");

        // theta
        // clock.tick("gerar_theta");
        A = idxA.size();
        B = idxB.size();
        theta(it) = R::rbeta(1 + A, 1 + B);
        // clock.tock("gerar_theta");

        // phi
        // clock.tick("subset_x");
        XA = x.rows(idxA);
        XB = x.rows(idxB);
        // clock.tock("subset_x");
        // clock.tick("subset_log_y");
        yA = log_y.elem(idxA);
        yB = log_y.elem(idxB);
        // clock.tock("subset_log_y");

        // clock.tick("gerar_phi");
        phi_a(it) = gerar_phi(A, yA, XA, aux_beta_a);
        phi_b(it) = gerar_phi(B, yB, XB, aux_beta_b);
        // clock.tock("gerar_phi");

        // clock.tick("gerar_beta");
        beta_a.col(it) = gerar_beta(yA, XA, phi_a(it));
        beta_b.col(it) = gerar_beta(yB, XB, phi_b(it));
        // clock.tock("gerar_beta");
    }
    // clock.tock("algoritmo_completo");
    // clock.stop("run_times");
    return Rcpp::List::create(
        Rcpp::_["betaA"] = beta_a,
        Rcpp::_["betaB"] = beta_b,
        Rcpp::_["phiA"] = phi_a,
        Rcpp::_["phiB"] = phi_b,
        Rcpp::_["theta"] = theta);
}
