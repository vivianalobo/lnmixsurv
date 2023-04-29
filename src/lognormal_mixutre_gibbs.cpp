// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// [[Rcpp::plugins("cpp11")]]

#include <RcppArmadillo.h>
#include <mvnorm.h>
#include <truncnorm.h>

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

// Realiza n amostras de uma dsitribuicao categorica.
// prob deve ser uma matriz n x k, onde n é o numero de observacoes
// e k o numero de categorias.
arma::colvec rcat(int n, arma::mat prob)
{
    int k = prob.n_cols;
    arma::colvec unif(n, arma::fill::randu);
    arma::colvec ret(n);
    arma::mat cumsum_prob = arma::cumsum(prob, 1);
    for (int obs = 0; obs < n; obs++)
    {
        for (int i = 0; i < k; i++)
        {
            if (unif(obs) < cumsum_prob(obs, i))
            {
                ret(obs) = i;
                break;
            }
        }
    }
    return ret;
}


// rtn1 mas com alguns argumentos especificos vetorizados.
arma::colvec rtn1_vec(arma::uword n, arma::colvec mean, double sd, arma::colvec low, double high)
{
    arma::colvec ret(n);
    for (arma::uword i = 0; i < n; i++)
    {
        ret(i) = r_truncnorm(mean(i), sd, low(i), high);
    }
    return ret;
}

// Random number generator for a dirichlet distribution.
// source: https://www.mjdenny.com/blog.html
arma::mat rdirichlet(int num_samples, arma::vec alpha_m)
{
    int distribution_size = alpha_m.n_elem;
    // each row will be a draw from a Dirichlet
    arma::mat distribution = arma::zeros(num_samples, distribution_size);

    for (int i = 0; i < num_samples; ++i)
    {
        double sum_term = 0;
        // loop through the distribution and draw Gamma variables
        for (int j = 0; j < distribution_size; ++j)
        {
            double cur = R::rgamma(alpha_m[j], 1.0);
            distribution(i, j) = cur;
            sum_term += cur;
        }
        // now normalize
        for (int j = 0; j < distribution_size; ++j)
        {
            distribution(i, j) = distribution(i, j) / sum_term;
        }
    }
    return (distribution);
}

arma::mat calcular_prob(arma::colvec log_y, arma::rowvec theta, arma::mat x, arma::mat beta, arma::rowvec sd)
{
    // n = numero de observacoes
    // p = numero de covariaveis
    // k = numero de componentes
    // n x p * p x k = n x k
    int n = log_y.size();
    int k = theta.size();
    arma::mat mean = x * beta;
    arma::mat prob(n, k);
    for (int i = 0; i < k; i++)
    {
        prob.col(i) = theta(i) * dnorm_vec(log_y, mean.col(i), sd(i), false);
    }
    arma::colvec sum_prob = arma::sum(prob, 1);    
    return prob / arma::repmat(sum_prob, 1, k);
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
arma::field<arma::cube> lognormal_mixture_gibbs(
    arma::mat x, arma::colvec y, arma::colvec delta,
    arma::uword iter = 1000, double valor_inicial_beta = 0,
    int numero_componentes_mistura = 2)
{

    int numero_observacoes = y.size();
    int numero_covariaveis = x.n_cols;

    arma::field<arma::cube> ret(3);
    arma::colvec log_y = log(y);
    arma::colvec c = arma::colvec(log_y);
    arma::uvec cens = (delta == 0);
    arma::cube theta(1, numero_componentes_mistura, iter);
    arma::cube phi(1, numero_componentes_mistura, iter);
    arma::cube beta(numero_covariaveis, numero_componentes_mistura, iter);

    // valores iniciais
    theta.slice(0) = rdirichlet(1, arma::vec(numero_componentes_mistura, arma::fill::ones));
    for (int i = 0; i < numero_componentes_mistura; i++) {
        phi.slice(0).col(i) = R::runif(0, 1);
    }

    beta.slice(0) = arma::mat(numero_covariaveis, numero_componentes_mistura, arma::fill::value(valor_inicial_beta));

    // variaveis auxiliares usadas no for
    arma::rowvec sd;
    arma::mat prob(numero_observacoes, numero_componentes_mistura);
    arma::colvec I(numero_observacoes);
    arma::field<arma::uvec> idx(numero_componentes_mistura);
    arma::colvec idx_sizes(numero_componentes_mistura);
    arma::mat beta_corrente(numero_covariaveis, numero_componentes_mistura);

    for (arma::uword it = 1; it <= iter - 1; it++) {

        beta_corrente = beta.slice(it - 1);

        sd = 1 / sqrt(phi.slice(it - 1));

        // probability
        prob = calcular_prob(log_y, theta.slice(it - 1), x, beta_corrente, sd);
        
        // mixture
        I = rcat(numero_observacoes, prob);


        for (int k = 0; k < numero_componentes_mistura; k++) {
            realizar_augmentation(log_y, c, x, cens, beta_corrente.col(k), sd(k), I, k);
            idx(k) = arma::find(I == k);
            idx_sizes(k) = idx(k).size();
        }
        
        // theta

        theta.slice(it) = rdirichlet(1, 1 + idx_sizes);

        // phi
        for (int k = 0; k < numero_componentes_mistura; k++) {
            phi.slice(it).col(k) = gerar_phi(idx_sizes(k), log_y.elem(idx(k)), x.rows(idx(k)), beta_corrente.col(k));
        }

        for (int k = 0; k < numero_componentes_mistura; k++) {
            beta.slice(it).col(k) = gerar_beta(log_y.elem(idx(k)), x.rows(idx(k)), phi.slice(it)(0, k));
        }
        
    }
    ret(0) = beta;
    ret(1) = phi;
    ret(2) = theta;
    return ret;
}


// [[Rcpp::export]]
arma::field<arma::field<arma::cube>> sequential_lognormal_mixture_gibbs(
    arma::mat x, arma::colvec y, arma::colvec delta,
    int iter, int chains, double valor_inicial_beta = 0, int numero_componentes = 2)
{
    arma::field<arma::field<arma::cube>> ret(chains);

    for (int i = 0; i < chains; i++)
    {
        try
        {
            ret(i) = lognormal_mixture_gibbs(x, y, delta, iter, valor_inicial_beta, numero_componentes);
        }
        catch (...)
        {
            Rcpp::warning("One or more chains were not generated because of some error.");
        }
    }
    return ret;
}
