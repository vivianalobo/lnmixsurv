// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// [[Rcpp::plugins("cpp11")]]

#include <RcppArmadillo.h>
#include <mvnorm.h>
#include "RcppTN.h"
// #include <progress.hpp>
// #include <progress_bar.hpp>
#include <omp.h> 

// #include <RcppClock.h>

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
    arma::colvec unif(n, arma::fill::randu);
    arma::uvec ret = unif < prob;
    return arma::conv_to<arma::vec>::from(ret);
}

// rtn1 mas com alguns argumentos especificos vetorizados.
arma::colvec rtn1_vec(arma::uword n, arma::colvec mean, double sd, arma::colvec low, double high)
{
    arma::colvec ret(n);
    for (arma::uword i = 0; i < n; i++)
    {
        ret(i) = RcppTN::rtn1(mean(i), sd, low(i), high);
    }
    return ret;
}

arma::colvec calcular_prob(arma::colvec log_y, double theta, arma::mat x, arma::mat beta, double sd_a, double sd_b)
{
    // n x p * p x 2 = n x 2
    arma::mat mean = x * beta;
    arma::colvec mean_a = mean.col(0);
    arma::colvec mean_b = mean.col(1);
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
//'
//' @noRd
void realizar_augmentation(arma::colvec& log_y, arma::colvec c, arma::mat x, arma::uvec cens, arma::colvec t_beta, double sd, arma::colvec indicadora_grupo, int grupo)
{
    arma::uvec g_z = arma::find(indicadora_grupo == grupo && cens == 1);

    // augmentation
    arma::colvec g_z_mean = x.rows(g_z) * t_beta;

    log_y.elem(g_z) = rtn1_vec(g_z.size(), g_z_mean, sd, c.elem(g_z), R_PosInf);
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
//' @noRd
arma::field<arma::mat> lognormal_mixture_gibbs_cpp(
    arma::mat x, arma::colvec y, arma::colvec delta,
    arma::uword iter = 1000, double valor_inicial_beta = 0)
{

    int numero_observacoes = y.size();
    int numero_covariaveis = x.n_cols;
    int numero_componentes_mistura = 2;

    arma::field<arma::mat> ret(5);
    arma::colvec log_y = log(y);
    arma::uvec cens = (delta == 0);
    arma::colvec theta(iter);
    arma::colvec phi_a(iter);
    arma::colvec phi_b(iter);
    arma::cube beta(numero_covariaveis, numero_componentes_mistura, iter);
    arma::colvec c = arma::colvec(log_y);

    // valores iniciais
    theta(0) = R::runif(0, 1);
    phi_a(0) = R::runif(0, 1);
    phi_b(0) = R::runif(0, 1);

    beta.slice(0) = arma::mat(numero_covariaveis, numero_componentes_mistura, arma::fill::value(valor_inicial_beta));

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

    arma::mat beta_corrente(numero_covariaveis, numero_componentes_mistura);

    // profiling
    // Rcpp::Clock clock;

    // Progress pro(iter, true);
    // clock.tick("algoritmo_completo");
    for (arma::uword it = 1; it <= iter - 1; it++)
    {
        // clock.tick("check_abort");
        // if (Progress::check_abort())
        // {
        //     // clock.stop("run_times");
        //     return Rcpp::List::create(
        //         Rcpp::_["beta"] = beta,
        //         Rcpp::_["phi_a"] = phi_a,
        //         Rcpp::_["phi_b"] = phi_b,
        //         Rcpp::_["theta"] = theta);
        // }

        // clock.tock("check_abort");
        // clock.tick("increment_pb");
        // pro.increment();
        // clock.tock("increment_pb");

        // clock.tick("criar_aux_beta");
        beta_corrente = beta.slice(it - 1);
        aux_beta_a = beta_corrente.col(0);
        aux_beta_b = beta_corrente.col(1);
        // clock.tock("criar_aux_beta");

        // clock.tick("criar_sd");
        sd_a = 1 / sqrt(phi_a(it - 1));
        sd_b = 1 / sqrt(phi_b(it - 1));
        // clock.tock("criar_sd");

        // probability
        // clock.tick("calcular_prob");
        // Rcpp::Rcout << "Antes prob" << std::endl;
        prob = calcular_prob(log_y, theta(it - 1), x, beta_corrente, sd_a, sd_b);
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
        // Rcpp::Rcout << "Antes aug" << std::endl;
        realizar_augmentation(log_y, c, x, cens, aux_beta_a, sd_a, I, 1);
        realizar_augmentation(log_y, c, x, cens, aux_beta_b, sd_b, I, 0);
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
        // Rcpp::Rcout << "Antes phi" << std::endl;
        phi_a(it) = gerar_phi(A, yA, XA, aux_beta_a);
        phi_b(it) = gerar_phi(B, yB, XB, aux_beta_b);
        // clock.tock("gerar_phi");

        // clock.tick("gerar_beta");
        // Rcpp::Rcout << "Antes beta" << std::endl;
        beta.slice(it).col(0) = gerar_beta(yA, XA, phi_a(it));
        beta.slice(it).col(1) = gerar_beta(yB, XB, phi_b(it));
        // clock.tock("gerar_beta");
    }
    // clock.tock("algoritmo_completo");
    // clock.stop("run_times");
    ret(0) = arma::mat(beta.col(0)).t();
    ret(1) = arma::mat(beta.col(1)).t();
    ret(2) = phi_a;
    ret(3) = phi_b,
        ret(4) = theta;
    return ret;
}

// [[Rcpp::export]]
arma::field<arma::cube> parallel_lognormal_mixture_gibbs_cpp(
    arma::mat x, arma::colvec y, arma::colvec delta,
    int iter, int chains, int cores = 1, double valor_inicial_beta = 0)
{

    arma::field<arma::mat> current_post(5);
    arma::field<arma::cube> ret(5);
    arma::cube beta_a(iter, chains, x.n_cols);
    arma::cube beta_b(iter, chains, x.n_cols);
    arma::cube phi_a(iter, chains, 1);
    arma::cube phi_b(iter, chains, 1);
    arma::cube theta(iter, chains, 1);
    omp_set_dynamic(0);

#pragma omp parallel for private(current_post) shared(y, x, delta, valor_inicial_beta, ret, iter) num_threads(cores)
    for (int i = 0; i < chains; i++) {

        bool exception_caught = true;
        try
        {
            current_post = lognormal_mixture_gibbs_cpp(x, y, delta, iter, valor_inicial_beta);
            exception_caught = false;
        }
        catch (...)
        {
            std::cout << "Erro na cadeia " << i << "\n";
        }
        if (!exception_caught) {
            beta_a.col(i) = current_post(0);
            beta_b.col(i) = current_post(1);
            phi_a.col(i) = current_post(2);
            phi_b.col(i) = current_post(3);
            theta.col(i) = current_post(4);
        }
    }
    ret(0) = beta_a;
    ret(1) = beta_b;
    ret(2) = phi_a;
    ret(3) = phi_b;
    ret(4) = theta;
    return ret;
}

// TODO: Não deveria ter essa repeticao de codigo com relacao ao parallel,
// mas o parallel com 1 core está demorando absurdos e nao to conseguindo
// resolver. Talvez depois tester com C++ puro para ver se o problema 
// aparece lá tb ou só no R
// [[Rcpp::export]]
arma::field<arma::cube> sequential_lognormal_mixture_gibbs_cpp(
    arma::mat x, arma::colvec y, arma::colvec delta,
    int iter, int chains, double valor_inicial_beta = 0)
{

    arma::field<arma::mat> current_post(5);
    arma::field<arma::cube> ret(5);
    arma::cube beta_a(iter, chains, x.n_cols);
    arma::cube beta_b(iter, chains, x.n_cols);
    arma::cube phi_a(iter, chains, 1);
    arma::cube phi_b(iter, chains, 1);
    arma::cube theta(iter, chains, 1);

    for (int i = 0; i < chains; i++) {

        bool exception_caught = true;
        try
        {
            current_post = lognormal_mixture_gibbs_cpp(x, y, delta, iter, valor_inicial_beta);
            exception_caught = false;
        }
        catch (...)
        {
            std::cout << "Erro na cadeia " << i << "\n";
        }
        if (!exception_caught) {
            beta_a.col(i) = current_post(0);
            beta_b.col(i) = current_post(1);
            phi_a.col(i) = current_post(2);
            phi_b.col(i) = current_post(3);
            theta.col(i) = current_post(4);
        }
    }
    ret(0) = beta_a;
    ret(1) = beta_b;
    ret(2) = phi_a;
    ret(3) = phi_b;
    ret(4) = theta;
    return ret;
}