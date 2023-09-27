// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#include <omp.h>
#include <RcppArmadillo.h>
#include <RcppGSL.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>

// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;

// ------ RNG Framework ------
// Function to initialize the GSL random number generator
void initializeRNG(const long long int& seed, gsl_rng* rng_device) {
  gsl_rng_set(rng_device, seed);
}

void setSeed(const long long int& seed, gsl_rng* rng_device) {
  initializeRNG(seed, rng_device);
}

double runif_0_1(gsl_rng* rng_device) {
  return gsl_rng_uniform(rng_device);
}

double rnorm_(const double& mu, const double& sd, gsl_rng* rng_device) {
  return gsl_ran_gaussian(rng_device, sd) + mu;
}

double rgamma_(const double& alpha, const double& beta, gsl_rng* rng_device) {
  return gsl_ran_gamma(rng_device, alpha, 1.0 / beta);
}

// Sample one value (k-dimensional) from a 
// Dirichlet(alpha_1, alpha_2, ..., alpha_k)
arma::vec rdirichlet(const arma::vec& alpha, gsl_rng* rng_device) {
  int K = alpha.n_elem;
  arma::vec sample(K);
  
  for (int k = 0; k < K; ++k) {
    sample(k) = rgamma_(alpha(k), 1.0, rng_device);
  }
  
  sample /= arma::sum(sample);
  return sample;
}

arma::vec rmvnorm(const arma::vec& mean, const arma::mat& covariance, gsl_rng* rng_device) {
  int numDims = mean.n_elem;
  arma::vec sample(numDims);
  
  arma::mat L = arma::chol(covariance, "lower");
  
  arma::vec Z(numDims);
  
  for (int j = 0; j < numDims; j++) {
    Z(j) = rnorm_(0, 1, rng_device);
  }
  
  sample = mean + L * Z;
  return sample;
}

/* AUXILIARY FUNCTIONS */
// Creates a sequence from start to end with 1 step
arma::ivec seq(const int& start, const int& end) {
  arma::vec out_vec = arma::linspace<arma::vec>(start, end, end - start + 1);
  return arma::conv_to<arma::ivec>::from(out_vec);
}

// Function for replicating a numeric value K times.
arma::vec repl(const double& x, const int& times) {
  return arma::ones<arma::vec>(times) * x;
}

// Sample a random object from a given vector
// for now, just sample numeric objects (because of c++ class definition)
// and just one object.
int numeric_sample(const arma::ivec& groups,
                   const arma::vec& probs, gsl_rng* rng_device) {
  double u = runif_0_1(rng_device);
  double cumulativeProb = 0.0;
  int n = probs.n_elem;
  for (int i = 0; i < n; ++i) {
    cumulativeProb += probs(i);
    
    if (u <= cumulativeProb) {
      return groups(i);
    }
  }
  
  // This point should never be reached and it's here just for compiling issues
  return 0;
}

// Function used to sample the latent groups for each observation. This one is
// the one that causes the code to run SO slow. We can further investigate what
// can be done in here
arma::ivec sample_groups(const int& G, const arma::vec& y, const arma::vec& eta, 
                         const arma::vec& phi, const arma::mat& beta,
                         const arma::mat& X, gsl_rng* rng_device) {
  arma::ivec vec_grupos(y.size());
  
  arma::vec probs(G);
  
  arma::mat mean = X * beta.t();
  arma::vec sd = 1.0 / sqrt(phi);
  double denom;
  int n = y.n_elem;
  for (int i = 0; i < n; i++) {
    denom = 0;
    
    for (int g = 0; g < G; g++) {
      probs(g) = eta(g) * R::dnorm(y(i), mean(i, g), sd(g), false);
      denom += probs(g); 
    }
    
    probs = (denom == 0) * (repl(1.0 / G, G)) +
      (denom != 0) * (probs / denom);
    
    vec_grupos(i) = numeric_sample(seq(0, G - 1), probs, rng_device);
  }
  
  return vec_grupos;
}

// Function used to simulate survival time for censored observations.
// Here, delta is a vector such that delta(i) == 0 implies that the 
// i-th observation is censored. Otherwise, the i-th observation is a
// failure time.

// This function can be vectorized with a little of work, which can
// improve performance. Take a look at
// https://cran.r-project.org/web/packages/RcppTN/RcppTN.pdf
// for sampling multiples numbers at the same time.
arma::vec augment(int G, const arma::vec& y, const arma::ivec& groups,
                  const arma::ivec& delta, const arma::vec& phi, 
                  const arma::mat& beta, const arma::mat& X, gsl_rng* rng_device) {
  
  int n = X.n_rows;
  
  arma::vec out(n);
  int g;
  double out_i;
  
  for (int i = 0; i < n; i++) {
    if (delta(i) == 1) {
      out(i) = y(i);
    }
    
    else {
      g = groups(i);
      
      out_i = y(i);
      int count = 0;
      while(out_i <= y(i)) {
        out_i = rnorm_(arma::as_scalar(X.row(i) * beta.row(g).t()),
                       sqrt(1.0 / phi(g)), rng_device);
        
        // break if it's going to run forever
        // should never happen
        if(count > 100000) {
          out_i = y(i) + 0.1;
          break;
        }
        count ++;
      }
      
      out(i) = out_i;
    }
  }
  
  return out;
}

// Create a table for each numeric element in the vector groups.
// For now, we are just going to use to evaluate how much observations
// we have at each group, given a vector of latent groups.
arma::ivec groups_table(const int& G, const arma::ivec& groups) {
  arma::ivec out(G);
  arma::ivec index;
  for (int g = 0; g < G; g++) {
    index = groups(arma::find(groups == g));
    
    out(g) = index.n_rows;
  }
  
  return(out);
}

arma::mat makeSymmetric(const arma::mat X) {
  arma::mat out(X.n_rows, X.n_cols);
  int rows = X.n_rows;
  int cols = X.n_cols;
  for(int r = 0; r < rows; r++) {
    for(int c = 0; c < cols; c++) {
      out(r, c) = X(r, c);
      out(c, r) = X(r, c);
    }
  }
  
  return out;
}

arma::mat lognormal_mixture_gibbs_implementation(int Niter, int em_iter, int G, 
                                                 arma::vec exp_y, arma::ivec delta, 
                                                 arma::mat X, double a, 
                                                 long long int starting_seed,
                                                 bool show_output, int chain_num) {
  
  gsl_rng* global_rng = gsl_rng_alloc(gsl_rng_default);
  
  // setting global seed to start the sampler
  setSeed(starting_seed, global_rng);
  
  // add verifications for robustiness. Skipping for the sake of simplicity.
  
  // Calculating number of columns of the output matrix:
  // Each group has p (#cols X) covariates, 1 mixture component and
  // 1 precision. This implies:
  int p = X.n_cols;
  int nColsOutput = (p + 2) * G;
  int N = X.n_rows;
  
  arma::vec y = log(exp_y);
  
  // The output matrix should have Niter rows (1 row for each iteration) and
  // nColsOutput columns (1 column for each element).
  arma::mat out(Niter, nColsOutput);
  
  // The order of filling the output matrix matters a lot, since we can
  // make label switching accidentally. Latter this is going to be defined
  // so we can always fill the matrix in the correct order (by columns, always).
  
  // Starting empty objects for EM
  arma::vec eta_em(G);
  arma::mat Xt = X.t();
  arma::vec phi_em(G);
  arma::vec y_aug(N);
  arma::mat beta_em(G, p);
  arma::mat tau(N, G);
  arma::sp_mat Wg(N, N);
  arma::ivec n_groups(G);
  arma::mat means(N, G);
  arma::vec sd(G);
  arma::vec linearComb(N);
  
  double sumtau;
  
  // starting EM algorithm to find values close to the MLE
  for (int iter = 0; iter < em_iter; iter++) {
    
    if ((iter % 50 == 0) && show_output) {
      Rcout << "(Chain " << chain_num << ") EM Iter: " << iter << "/" << em_iter << "\n";
    }
    
    // Initializing values
    if(iter == 0) {
      eta_em = rdirichlet(repl(1, G), global_rng);
      
      for (int g = 0; g < G; g++) {
        phi_em(g) = rgamma_(2, 8, global_rng);
        beta_em.row(g) = rmvnorm(repl(0, p),
                    arma::diagmat(repl(7, p)),
                    global_rng).t();
      }
    }
    
    // E-step
    means = X * beta_em.t();
    sd = sqrt(1.0 / phi_em);
    
    for(int r = 0; r < N; r++) {
      for(int g = 0; g < G; g++) {
        tau(r, g) = eta_em(g) * R::dnorm(y(r), 
            means(r, g), sd(g), false);
      }
      
      if(arma::sum(tau.row(r)) == 0) {  // to avoid numerical problems
        tau.row(r) = repl(1.0/G, G).t();
      } else {
        tau.row(r) = tau.row(r)/sum(tau.row(r));
      }
    }
    
    // M-step
    for(int g = 0; g < G; g++) {
      Wg = arma::diagmat(tau.col(g));
      sumtau = arma::sum(tau.col(g));
      
      eta_em(g) = sumtau/N;
      
      if(arma::det(X.t() * Wg * X) != 0) {
        beta_em.row(g) = (arma::inv(Xt * Wg * X,
                          arma::inv_opts::allow_approx) * Xt * Wg * y).t();
      }
      
      linearComb = y - X * beta_em.row(g).t();
      
      phi_em(g) = sumtau/arma::as_scalar(linearComb.t() * Wg * linearComb);
    }
  }
  
  // Starting other new values for MCMC algorithms
  arma::vec eta(G);
  arma::vec phi(G);
  arma::mat beta(G, p);
  arma::ivec groups(y.n_rows);
  arma::vec e0new_vec(2);
  
  arma::mat Xg;
  arma::mat Xgt;
  arma::vec yg;
  
  arma::uvec indexg;
  arma::mat Sg(p, p);
  arma::vec mg(p);
  arma::vec log_eta_new(G);
  
  double e0;
  double e0_prop;
  double log_alpha;
  double b;
  double cte;
  double prop_aceite;
  int n_aceite;
  double count = 0;
  
  for (int iter = 0; iter < Niter; iter++) {
    // Starting empty objects for Gibbs Sampler
    if (iter == 0) { 
      // we are going to start the values using the
      // previous EM iteration
      eta = eta_em;
      phi = phi_em;
      beta = beta_em;
      
      // sampling value for e0
      e0 = rgamma_(1, 1, global_rng);
      
      // defining values for sintonizing the variance
      // of e0 proposal
      cte = 1;
      n_aceite = 0;
      
      // Sampling classes for the observations
      groups = sample_groups(G, y, eta, phi, beta, X, global_rng);
    }
    
    // Data augmentation
    y_aug = augment(G, y, groups, delta, phi, beta, X, global_rng);
    
    // Sampling classes for the observations
    groups = sample_groups(G, y_aug, eta, phi, beta, X, global_rng);
    
    // Computing number of observations allocated at each class
    n_groups = groups_table(G, groups);
    
    // ensuring that every class have, at least, 2 observations
    for(int g = 0; g < G; g++) {
      if(n_groups(g) == 0) {
        for(int m = 0; m < 2; m++) {
          int idx = numeric_sample(seq(0, X.n_rows),
                                   repl(1.0 / X.n_rows, X.n_rows),
                                   global_rng);
          groups(idx) = g;
        }
        
        // recalculating the number of groups
        n_groups = groups_table(G, groups);
      }
    }
    
    // Sampling new eta
    eta = rdirichlet(arma::conv_to<arma::Col<double>>::from(n_groups) + e0, global_rng);
    
    // For each g, sample new phi[g] and beta[g, _]
    for (int g = 0; g < G; g++) {
      indexg = arma::find(groups == g);
      if(true) {
        Xg = X.rows(indexg);
        Xgt = Xg.t();
        yg = y_aug(indexg);
        
        linearComb = yg - Xg * beta.row(g).t();
        
        // sampling phi new
        // the priori used was Gamma(0.001, 0.001)
        
        phi(g) = rgamma_(n_groups(g) / 2.0 + 0.001, (1.0 / 2) *
          as_scalar(linearComb.t() * linearComb) + 0.001, global_rng);
        
        // sampling beta new
        // the priori used was MNV(vec 0, diag 1000)
        if(arma::det(phi(g) * Xgt * Xg + 
           arma::diagmat(repl(1.0 / 1000, p))) != 0) {
          Sg = arma::inv(phi(g) * Xgt * Xg + 
            arma::diagmat(repl(1.0 / 1000, p)));
          
          if(Sg.is_symmetric() == false) {
            Sg = makeSymmetric(Sg);
          }
          
          mg = Sg * (phi(g) * Xgt * yg);
          
          beta.row(g) = rmvnorm(mg, Sg, global_rng).t();
        }
      }
    }
    
    // atualizando o valor da constante de sintonização
    // cte a cada 200 iterações
    if ((iter % 200) == 0) {
      prop_aceite = n_aceite/200.0;
      
      if (prop_aceite > 0.25) {
        cte = cte/((2*sqrt(3)-2)/(1+exp(0.04*count)) + 1);
      } else if (prop_aceite < 0.17) {
        cte = cte*((2*sqrt(3)-2)/(1+exp(0.04*count)) + 1);
      }
      // ((2*sqrt(3)-2)/(1+exp(0.04*count)) + 1) é um valor
      // de correção que converge para 1 quando count -> Inf e,
      // quando count = 0, o resultado é sqrt(3).
      
      // o valor 0.04 representa a velocidade da convergência
      // para 1 e sqrt(3), o valor em count = 0, foi definido
      // arbitrariamente.
      
      n_aceite = 0;
      count = count + 1.0;
    }
    
    b = cte*a;
    
    // updating the value of e0 (eta's dirichlet hyperparameter)
    e0_prop = rgamma_(b*e0, b, global_rng);
    
    log_eta_new = log(eta);
    
    log_alpha = arma::sum(log_eta_new)*(e0_prop - e0) +
      9.0 * log(e0_prop/e0) - 10.0*G*(e0_prop - e0) +
      b*(e0_prop - e0) + (b*e0_prop - 1.0)*log(e0) -
      (b*e0 - 1.0)*log(e0_prop);
    
    double u = runif_0_1(global_rng);
    e0 = (log(u) < log_alpha) * e0_prop +
      (log(u) >= log_alpha) * e0;
    
    n_aceite += (log(u) < log_alpha);
    
    // filling the ith iteration row of the output matrix
    // the order of filling will always be the following:
    
    // First Mixture: proportion, betas, phi
    // Second Mixture: proportion, betas, phi
    // ...
    // Last Mixture: proportion, betas, phi
    
    // arma::uvec sorteta = arma::sort_index(eta, "descend");
    // beta = beta.rows(sorteta);
    // phi = phi.rows(sorteta);
    // eta = eta.rows(sorteta);
    
    arma::rowvec newRow = arma::join_rows(beta.row(0),
                                          phi.row(0),
                                          eta.row(0));
    for (int g = 1; g < G; g++) {
      newRow = arma::join_rows(newRow, beta.row(g),
                               phi.row(g),
                               eta.row(g));
    }
    
    out.row(iter) = newRow;
    
    if((iter % 500 == 0) && show_output) {
      Rcout << "(Chain " << chain_num << ") MCMC Iter: " << iter << "/" << Niter << "\n";
    }
  }
  
  if(show_output) {
    Rcout << "Chain " << chain_num << " finished sampling." << "\n";
  }
  return out;
}

// Function to call lognormal_mixture_gibbs_implementation with parallellization
// [[Rcpp::export]]
arma::cube lognormal_mixture_gibbs(int Niter, int em_iter, int G, arma::vec exp_y, arma::ivec delta, arma::mat X,
                                   double a, arma::Col<long long int> starting_seed, bool show_output, int n_cores, int n_chains,
                                   bool force_num_cores) {
  arma::cube out(Niter, (X.n_cols + 2) * G, n_chains);
  
  if(n_cores == 1) {
    for(int chain = 0; chain < n_chains; chain ++) {
      out.slice(chain) = lognormal_mixture_gibbs_implementation(Niter, em_iter, G, exp_y, delta, X, a, starting_seed(chain), show_output,
                chain + 1);
    }
    
    return out;
  }
  
  if(force_num_cores) {
    omp_set_dynamic(0); // related to https://stackoverflow.com/questions/11095309/openmp-set-num-threads-is-not-working
  }
  
  omp_set_num_threads(n_cores);
  
  int chain;
  
  #pragma omp parallel for private(chain)
  for(int chain = 0; chain < n_chains; chain ++) {
    #pragma omp critical
    usleep(5000 * chain); // sleep to avoid racing conditions at the beginning
    out.slice(chain) = lognormal_mixture_gibbs_implementation(Niter, em_iter, G, exp_y, delta, X, a, 
              starting_seed(chain), show_output, chain + 1);
  }
  
  return out;
}