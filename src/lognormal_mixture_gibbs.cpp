// -*- mode: C++; c-indent-level: 2; c-basic-offset: 2; indent-tabs-mode: nil; -*-

#include <RcppArmadillo.h>
#include <RcppGSL.h>
#include <RcppParallel.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>

using namespace Rcpp;

// Importing the RcppParallelLibs Function from RcppParallel Package to NAMESPACE
//' @importFrom RcppParallel RcppParallelLibs
 
 // ------ RNG Framework ------
 
 // Function to initialize the GSL random number generator
 void initializeRNG(const long long int& seed, gsl_rng* rng_device) {
   gsl_rng_set(rng_device, seed);
 }

// Function used to set a seed
void setSeed(const long long int& seed, gsl_rng* rng_device) {
  initializeRNG(seed, rng_device);
}

// Squares a double
double square(const double& x) {
  return x * x;
}

// Generates a random observation from Uniform(0, 1) 
double runif_0_1(gsl_rng* rng_device) {
  return gsl_rng_uniform(rng_device);
}

// Generates a random observation from Normal(mu, sd^2)
double rnorm_(const double& mu, const double& sd, gsl_rng* rng_device) {
  return gsl_ran_gaussian(rng_device, sd) + mu;
}

// Generates a random observation from Gamma(alpha, beta), with mean alpha/beta
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

// Generates a random observation from a MultivariateNormal(mean, covariance)
arma::vec rmvnorm(const arma::vec& mean, const arma::mat& covariance, gsl_rng* rng_device) {
  int numDims = mean.n_elem;
  arma::vec sample(numDims);
  
  arma::mat L = arma::chol(covariance, "lower");
  
  arma::vec Z(numDims);
  
  for (int j = 0; j < numDims; j++) {
    Z(j) = rnorm_(0.0, 1.0, rng_device);
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
// Note: it just samples numeric objects (because of c++ class definition) and just one object per time.
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

/* NON SPARSE ALTERNATIVES */

// Function used to sample the latent groups for each observation.
arma::ivec sample_groups(const int& G, const arma::vec& y, const arma::vec& eta, 
                         const arma::vec& phi, const arma::mat& beta,
                         const arma::mat& X, gsl_rng* rng_device,
                         const arma::ivec& groups_old) {
  
  // Initializing variables used for sampling groups
  int n = X.n_rows;
  int p = X.n_cols;
  arma::ivec vec_groups = groups_old;
  arma::ivec sequence = seq(0, G - 1);
  arma::uvec indexg;
  arma::mat Vg0 = arma::diagmat(repl(30.0, p));
  arma::mat identity_p = arma::eye(p, p);
  arma::mat Vg0_inv = arma::solve(Vg0, identity_p, arma::solve_opts::allow_ugly);
  arma::vec xi;
  arma::rowvec xit;
  arma::mat Xg;
  arma::mat Xgt;
  arma::mat Xgt_Xg;
  arma::vec yg;
  arma::mat S_inv(p, p);
  arma::mat S(p, p);
  arma::vec Xgt_yg;
  arma::vec Mg;
  arma::mat denom_mat(n, G);
  arma::vec probsEqual(G);
  probsEqual.fill(1.0 / G);
  
  double m;
  double yi;
  double a;
  double sigma2;
  double denom;
  
  for (int g = 0; g < G; g++) {
    indexg = arma::find(vec_groups == g); // find indexes of which observations belongs to group g
    
    if(indexg.empty()) {
      continue; // continue to next group because if there's no one at group g, a lot of computational problem is going to happen next
    }
    
    Xg = X.rows(indexg);
    Xgt = Xg.t();
    yg = y(indexg);
    Xgt_yg = Xgt * yg;
    S = phi(g) * Xg.t() * Xg + Vg0_inv;
    
    if(arma::det(S) == 0) {
      S_inv = Vg0;
    } else {
      S_inv = arma::solve(S, identity_p, 
                          arma::solve_opts::allow_ugly);
    }
    
    Mg = X * S_inv * Xgt_yg;
    
    for (int i = 0; i < n; i++) {
      xit = X.row(i);
      yi = y(i);
      m = 1.0 / arma::as_scalar(xit * S_inv * xit.t());
      a = (m - phi(g)) * (vec_groups(i) == g) +
        m * (vec_groups(i) != g);
      
      sigma2 = (a + phi(g))/(a * phi(g));
      
      denom_mat(i, g) = 
        eta(g) * (R::dnorm(yi,
                  square(phi(g)) * sigma2 * (Mg(i) - yi / (a + phi(g))),
                  sqrt(sigma2), 
                  false) * (vec_groups(i) == g) +
                    R::dnorm(yi,
                             square(phi(g)) * sigma2 * ((a + phi(g)) * Mg(i)/(a + phi(g) + 1.0)),
                             sqrt(sigma2),
                             false) * (vec_groups(i) != g));
      
      if(g == (G - 1)) { // if it's the last group
        denom = arma::sum(denom_mat.row(i));
        
        if(denom > 0) {
          vec_groups(i) = numeric_sample(sequence, denom_mat.row(i).t() / denom, 
                     rng_device);
        } else {
          vec_groups(i) = numeric_sample(sequence, probsEqual, rng_device);
        }
      }
    }
  }
  
  return(vec_groups);
}

// Function used to sample random groups for each observation proportional to the eta parameter
arma::ivec sample_groups_start(const int& G, const arma::vec& y, 
                               const arma::vec& eta, gsl_rng* rng_device) {
  int n = y.n_rows;
  arma::ivec vec_groups(n);
  
  for (int i = 0; i < n; i++) {
    vec_groups(i) = numeric_sample(seq(0, G - 1), eta, rng_device);
  }
  
  return(vec_groups);
}

double augment_yi(const double& yi, const double& mean, const double& sd, gsl_rng* rng_device) {
  double out_i = yi; // value to augment
  int count = 0; // avoid infite loop
  
  while(out_i <= yi) { // ensure that the sampled observation is bigger than the censored one
    out_i = rnorm_(mean, sd, rng_device);
    
    // break if it seems like it's going to run forever
    // should never happen
    if(count > 10000) {
      out_i = yi + 0.01;
      break;
    }
    
    count ++;
  }
  
  return out_i;
}

// Function used to simulate survival time for censored observations.
// Here, delta is a vector such that delta(i) == 0 implies that the 
// i-th observation is censored. Otherwise, the i-th observation is a
// failure time.
arma::vec augment(const int& G, const arma::vec& y, const arma::ivec& groups,
                  const arma::ivec& delta, const arma::vec& sd, 
                  const arma::mat& beta, const arma::mat& X, gsl_rng* rng_device) {
  
  int n = X.n_rows;
  arma::vec out = y;
  arma::mat mean = X * beta.t(); // pre-compute the mean matrix
  arma::uvec censored_indexes = arma::find(delta == 0); // finding which observations are censored
  
  for (int i : censored_indexes) {
    out(i) = augment_yi(y(i), arma::as_scalar(mean(i, groups(i))), sd(groups(i)), rng_device);
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
  
  return out;
}

/*
 // Force matrix X to be symmetric. 
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
 */

/* Auxiliary functions for EM algorithm */

// Compute weights matrix
arma::mat compute_W(const arma::vec& y, const arma::mat& X, const arma::vec& eta, 
                    const arma::mat& beta, const arma::vec& sigma, 
                    const int& G) {
  
  int n = X.n_rows;
  double denom;
  arma::mat out(n, G);
  arma::mat mat_denom(n, G);
  arma::rowvec repl_vec = repl(1.0 / G, G).t();
  
  for(int g = 0; g < G; g++) {
    mat_denom.col(g) = eta(g) * arma::normpdf(y,
                  X * beta.row(g).t(),
                  repl(sigma(g), n));
  }
  
  for(int i = 0; i < n; i++) {
    denom = arma::sum(mat_denom.row(i));
    if(denom > 0) {
      out.row(i) = mat_denom.row(i) / denom;
    } else {
      out.row(i) = repl_vec;
    }
  }
  
  return(out);
}

// Function used to computed the expected value of a truncated normal distribution
double compute_expected_value_truncnorm(const double& alpha, const double& mean, const double& sigma) {
  double out;
  
  if (R::pnorm(alpha, 0.0, 1.0, true, false) < 1.0) {
    out = mean + sigma *
      (R::dnorm(alpha, 0.0, 1.0, false)/(1.0 - R::pnorm(alpha, 0.0, 1.0, true, false)));
  } else {
    out = mean + sigma *
      (R::dnorm(alpha, 0.0, 1.0, false)/(1.0 - 0.999));
  }
  
  return out;
}

// Create the latent variable z for censored observations
arma::vec augment_em(const arma::vec& y, const arma::ivec& delta,
                     const arma::mat& X, const arma::mat& beta,
                     const arma::vec& sigma, const arma::mat& W,
                     const int& G) {
  int n = X.n_rows;
  arma::vec out = y;
  arma::mat mean = X * beta.t();
  arma::mat alpha_mat(n, G);
  arma::uvec censored_indexes = arma::find(delta == 0); // finding which observations are censored
  
  for(int g = 0; g < G; g++) {
    alpha_mat.col(g) = (y - mean.col(g))/sigma(g);
  }
  
  for (int i : censored_indexes) {
    out(i) = 0.0;
    
    for (int g = 0; g < G; g++) {
      out(i) += W(i, g) * compute_expected_value_truncnorm(arma::as_scalar(alpha_mat(i, g)), arma::as_scalar(mean(i, g)), sigma(g));;
    }
  }
  
  return out;
}

// Function used to sample groups from W. It samples one group by row based on the max weight.
arma::ivec sample_groups_from_W(const arma::mat& W) {
  int N = W.n_rows;
  arma::vec out(N);
  
  for(int r = 0; r < N; r++) {
    out(r) = W.row(r).index_max();
  }
  
  return(arma::conv_to<arma::ivec>::from(out));
}

// Sample initial values for the EM parameters
void sample_initial_values_em(arma::vec& eta, arma::vec& phi, arma::mat& beta, arma::vec& sd, const int& G, const int& k, gsl_rng* rng_device) {
  eta = rdirichlet(repl(1.0, G), rng_device);
  
  for (int g = 0; g < G; g++) {
    phi(g) = rgamma_(0.5, 0.5, rng_device);
    
    for (int c = 0; c < k; c++) {
      beta(g, c) = rnorm_(0.0, 10.0, rng_device);
    }
  }
  
  sd = 1.0 / sqrt(phi);
}

// Update the matrix beta for the group g
void update_beta_g(const arma::vec& colg, const arma::mat& X, const int& g, const arma::vec& z, arma::mat& beta) {
  arma::sp_mat Wg;
  Wg = arma::diagmat(colg);
  
  if(arma::det(X.t() * Wg * X) != 0.0) {
    beta.row(g) = arma::solve(X.t() * Wg * X,
             X.t() * Wg * z,
             arma::solve_opts::allow_ugly).t();
  }
}

// Update the parameter phi(g)
void update_phi_g(const double& denom, const arma::uvec& censored_indexes, const arma::mat& X, const arma::vec& colg, const arma::vec& y, const arma::vec& z,
                  const arma::vec& sd, const arma::mat& beta, const arma::vec& var, const int& g, const int& n, arma::vec& phi, gsl_rng* rng_device) {
  double alpha = 0.0;
  double quant = arma::as_scalar(arma::square(z - (X * beta.row(g).t())).t() * colg);
  
  for(int i : censored_indexes) {
    alpha = (y(i) - arma::as_scalar(X.row(i) * beta.row(g).t())) / sd(g);
    
    if(R::pnorm(alpha, 0.0, 1.0, true, false) < 1.0) {
      quant += colg(i) * var(g) * (1.0 - 
        (- alpha * R::dnorm(alpha, 0.0, 1.0, false))/(1.0 - R::pnorm(alpha, 0.0, 1.0, true, false)) -
        square(R::dnorm(alpha, 0.0, 1.0, false)/(1.0 - R::pnorm(alpha, 0.0, 1.0, true, false))));
    } else {
      quant += colg(i) * var(g) * (1.0 - 
        (-alpha * R::dnorm(alpha, 0.0, 1.0, false))/(1.0 - 0.999) -
        square(R::dnorm(alpha, 0.0, 1.0, false)/(1.0 - 0.999)));
    }
  }
  
  // to avoid numerical problems
  if (quant == 0) {
    phi(g) = rgamma_(0.5, 0.5, rng_device); // resample phi
  }
  
  // to avoid numerical problems
  if(phi(g) > 1e5 || phi.has_nan()) {
    phi(g) = rgamma_(0.5, 0.5, rng_device); // resample phi
  }
}

// Update the model parameters with EM
void update_em_parameters(const int& n, const int& G, arma::vec& eta, arma::mat& beta, arma::vec& phi, const arma::mat& W, const arma::mat& X, 
                          const arma::vec& y, const arma::vec& z, const arma::ivec& delta, const arma::vec& sd, gsl_rng* rng_device) {
  arma::vec colg(n);
  arma::vec var;
  double quant, denom, alpha;
  arma::uvec censored_indexes = arma::find(delta == 0);
  var = arma::square(sd);
  
  for (int g = 0; g < G; g++) {
    colg = W.col(g);
    
    eta(g) = arma::sum(colg) / n; // updating eta(g)
    update_beta_g(colg, X, g, z, beta); // updating beta for the group g
    update_phi_g(arma::sum(colg), censored_indexes, X, colg, y, z, sd, beta, var, g, n, phi, rng_device);
  }
}

// Compute model's log-likelihood to select the EM initial values
double loglik_em(const arma::vec& eta, const arma::mat& beta, const arma::vec& phi, const arma::mat& W, const arma::vec& z, const int& G,
                 const arma::mat& X) {
  int N = X.n_rows;
  arma::vec sd = 1.0 / sqrt(phi);
  arma::mat mean = X * beta.t();
  double loglik = 0.0;
  
  for(int i = 0; i < N; i++) {
    for(int g = 0; g < G; g++) {
      loglik += W(i, g) * R::dnorm(z(i), arma::as_scalar(mean(i, g)),
                  sd(g), true);
    }
  }
  
  return loglik;
}

// EM for the lognormal mixture model.
arma::field<arma::mat> lognormal_mixture_em(const int& Niter, const int& G, const arma::vec& t, const arma::ivec& delta, const arma::mat& X,
                                            const bool& better_initial_values, const int& N_em,
                                            const int& Niter_em, const bool& internal, gsl_rng* rng_device) {
  
  int n = X.n_rows;
  int k = X.n_cols;
  
  // initializing objects used on EM algorithm
  arma::vec y = log(t);
  arma::vec eta(G);
  arma::vec phi(G);
  arma::vec sd(G);
  arma::vec z(n);
  arma::mat W(n, G);
  arma::mat beta(G, k);
  arma::mat out(Niter, G * k + (G * 2));
  arma::field<arma::mat> out_internal_true(6);
  arma::field<arma::mat> out_internal_false(2);
  
  arma::field<arma::mat> em_params(6);
  arma::field<arma::mat> best_em(6);
  for(int iter = 0; iter < Niter; iter++) {
    if(iter == 0) { // sample starting values
      if(better_initial_values) {
        for (int init = 0; init < N_em; init ++) {
          em_params = lognormal_mixture_em(Niter_em, G, t, delta, X, false, 0, 0, true, rng_device);
          
          if(init == 0) {
            best_em = em_params;
            Rcout << "Initial LogLik: " << arma::as_scalar(best_em(5)) << "\n";
          } else {
            if(arma::as_scalar(em_params(5)) > arma::as_scalar(best_em(5))) { // comparing logliks
              Rcout << "Previous maximum: " << arma::as_scalar(best_em(5)) << " | New maximum: " << arma::as_scalar(em_params(5))  << "\n";
              best_em = em_params;
            }
          }
        }
        
        eta = best_em(0);
        beta = best_em(1);
        phi = best_em(2);
        W = best_em(3);
      } else {
        sample_initial_values_em(eta, phi, beta, sd, G, k, rng_device);
        W = compute_W(y, X, eta, beta, sd, G);
      }
    } else {
      sd = 1.0 / sqrt(phi);
      z = augment_em(y, delta, X, beta, sd, W, G);
      W = compute_W(z, X, eta, beta, sd, G);
      update_em_parameters(n, G, eta, beta, phi, W, X, y, z, delta, sd, rng_device);
    }
    
    // Fill the out matrix
    arma::rowvec newRow = 
      arma::join_rows(eta.row(0), 
                      beta.row(0),
                      phi.row(0));
    
    for (int g = 1; g < G; g++) {
      newRow = 
        arma::join_rows(newRow, 
                        eta.row(g), 
                        beta.row(g),
                        phi.row(g));
    }
    
    out.row(iter) = newRow;
  }
  
  if(internal) {
    out_internal_true(0) = eta;
    out_internal_true(1) = beta;
    out_internal_true(2) = phi;
    out_internal_true(3) = W;
    out_internal_true(4) = augment_em(y, delta, X, beta, 1.0 / sqrt(phi), W, G);
    out_internal_true(5) = loglik_em(eta, beta, phi, W, out_internal_true(4), G, X);
    
    return out_internal_true;
  } else {
    out_internal_false(0) = out;
    out_internal_false(1) = loglik_em(eta, beta, phi, W, augment_em(y, delta, X, beta, 1.0 / sqrt(phi), W, G), G, X);
    
    return out_internal_false;
  }
  
  return out_internal_false; // should never be reached
}

// Setting parameter's values for the first Gibbs iteration
void first_iter_gibbs(const arma::field<arma::mat>& em_params, arma::vec& eta, arma::mat& beta, arma::vec& phi, const int& em_iter, const int& G, const arma::vec& y,
                      double& e0, arma::vec& sd, arma::ivec& groups, const arma::mat& X, const bool& use_W, const arma::ivec& delta, gsl_rng* rng_device) {
  arma::ivec groups_start(y.n_rows);
  int p = X.n_cols;
  if (em_iter != 0) {
    // we are going to start the values using the last EM iteration
    eta = em_params(0);
    beta = em_params(1);
    phi = em_params(2);
    sd = 1.0 / sqrt(phi);
    groups_start = sample_groups_from_W(em_params(3));
  } else {
    eta = rdirichlet(repl(1, G), rng_device);
    
    for (int g = 0; g < G; g++) {
      phi(g) = rgamma_(0.5, 0.5, rng_device);
      beta.row(g) = rmvnorm(repl(0.0, p),
               arma::diagmat(repl(15.0 * 15.0, p)),
               rng_device).t();
    }
    
    sd = 1.0 / sqrt(phi);
    // Sampling classes for the observations
    groups_start = sample_groups_start(G, y, eta, rng_device);
  }
  
  // sampling value for e0
  e0 = rgamma_(1.0, 1.0, rng_device);
  
  if(use_W == false) {
    groups = sample_groups(G, augment(G, y, groups_start, delta, sd, beta, X, rng_device), eta, phi, beta, X, rng_device, groups_start);
  } else {
    groups = groups_start;
  }
}

// Function to update groups
void update_groups_gibbs(const int& iter, const bool& use_W, const arma::field<arma::mat>& em_params, const int& G, const arma::vec& y_aug, 
                         const arma::vec& eta, const arma::mat& beta, const arma::vec& phi, const arma::mat& X, arma::ivec& groups, gsl_rng* rng_device) {
  if (iter > 0) {
    if(use_W) {
      groups = sample_groups_from_W(em_params(3));
    } else {
      groups = sample_groups(G, y_aug, eta, phi, beta, X, rng_device, groups);
    }
  }
}

// Avoiding groups with zero number of observations in it (causes numerical issues)
void avoid_group_with_zero_allocation(arma::ivec& n_groups, arma::ivec& groups, const int& G, const int& N, gsl_rng* rng_device) {
  int idx = 0;
  int m;
  
  for(int g = 0; g < G; g++) {
    if(n_groups(g) == 0) {
      m = 0;
      while(m < 5) {
        idx = numeric_sample(seq(0, N),
                             repl(1.0 / N, N),
                             rng_device);
        
        if(n_groups(groups(idx)) > 5) {
          groups(idx) = g;
          m += 1;
        } 
      }
      
      // recalculating the number of groups
      n_groups = groups_table(G, groups);
    }
  }
}

double update_phi_g_gibbs(const int& n_groups_g, const arma::vec& linearComb, gsl_rng* rng_device) {
  return rgamma_(n_groups_g / 2.0 + 0.01, (1.0 / 2.0) * arma::as_scalar(linearComb.t() * linearComb) + 0.01, rng_device);
}

arma::rowvec update_beta_g_gibbs(const double& phi_g, const arma::mat& Xg, const arma::mat& Xgt, const arma::vec& yg, gsl_rng* rng_device) {
  arma::rowvec out;
  arma::mat comb = phi_g * Xgt * Xg + arma::diagmat(repl(1.0 / 30.0, Xg.n_cols));
  arma::mat Sg;
  arma::vec mg;
  
  if(arma::det(comb) != 0) {
    Sg = arma::solve(comb,
                     arma::eye(Xg.n_cols, Xg.n_cols),
                     arma::solve_opts::allow_ugly);
    
    mg = arma::solve(comb,
                     phi_g * Xgt * yg,
                     arma::solve_opts::allow_ugly);
    
    out = rmvnorm(mg, Sg, rng_device).t();
  }
  
  return out;
}

// update all the Gibbs parameters
void update_gibbs_parameters(const int& G, const arma::mat& X, const arma::vec& y_aug, const arma::ivec& n_groups, const arma::ivec& groups, 
                             arma::vec& eta, arma::mat& beta, arma::vec& phi, double& e0, gsl_rng* rng_device) {
  
  arma::mat Xg;
  arma::mat Xgt;
  arma::vec yg;
  arma::vec linearComb;
  arma::uvec indexg;
  
  // updating eta
  eta = rdirichlet(arma::conv_to<arma::Col<double>>::from(n_groups) + e0, 
                   rng_device);
  
  // For each g, sample new phi[g] and beta[g, _]
  for (int g = 0; g < G; g++) {
    indexg = arma::find(groups == g);
    Xg = X.rows(indexg);
    Xgt = Xg.t();
    yg = y_aug.rows(indexg);
    linearComb = yg - Xg * beta.row(g).t();
    
    // updating phi(g)
    // the priori used was Gamma(0.01, 0.01)
    phi(g) = update_phi_g_gibbs(n_groups(g), linearComb, rng_device);
    
    // updating beta.row(g)
    // the priori used was MNV(vec 0, diag 1000)
    beta.row(g) = update_beta_g_gibbs(phi(g), Xg, Xgt, yg, rng_device);
  }
}

// update e0, hyperparameter of the eta's dirichlet prior
void update_e0(const int& G, const int& iter, double& prop_aceite, int& n_aceite, double& cte, const double& a, const arma::vec& eta, double& e0, double& count, gsl_rng* rng_device) {  
  arma::vec log_eta_new = log(eta);
  
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
  
  double b = cte * a;
  
  // updating the value of e0 (eta's dirichlet hyperparameter)
  double e0_prop = rgamma_(b*e0, b, rng_device);
  
  double log_alpha = arma::sum(log_eta_new) * (e0_prop - e0) + 9.0 * log(e0_prop/e0) -
    10.0*  G * (e0_prop - e0) + b * (e0_prop - e0) + (b * e0_prop - 1.0) *
    log(e0) - (b*e0 - 1.0) * log(e0_prop);
  
  double u = runif_0_1(rng_device);
  
  e0 = (log(u) < log_alpha) * e0_prop +
    (log(u) >= log_alpha) * e0;
  
  n_aceite += (log(u) < log_alpha);
}

// Internal implementation of the lognormal mixture model via Gibbs sampler
arma::mat lognormal_mixture_gibbs_implementation(const int& Niter, const int& em_iter, const int& G, 
                                                 const arma::vec& t, const arma::ivec& delta, 
                                                 const arma::mat& X, const double& a, 
                                                 long long int starting_seed,
                                                 const bool& show_output, const int& chain_num,
                                                 const bool& use_W, const bool& better_initial_values,
                                                 const int& Niter_em, const int& N_em) {
  
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
  
  arma::vec y = log(t);
  
  // The output matrix should have Niter rows (1 row for each iteration) and
  // nColsOutput columns (1 column for each element).
  arma::mat out(Niter, nColsOutput);
  
  // The order of filling the output matrix matters a lot, since we can
  // make label switching accidentally. Latter this is going to be defined
  // so we can always fill the matrix in the correct order (by columns, always).
  arma::mat Xt = X.t();
  arma::vec y_aug(N);
  arma::ivec n_groups(G);
  arma::mat means(N, G);
  arma::vec sd(G);
  
  // Starting other new values for MCMC algorithms
  arma::vec eta(G);
  arma::vec phi(G);
  arma::mat beta(G, p);
  arma::ivec groups(N);
  arma::vec log_eta_new(G);
  
  double e0;
  double prop_aceite;
  int n_aceite;
  double count = 0.0;
  double cte;

  arma::rowvec newRow;
  arma::field<arma::mat> em_params(6);
  
  if(em_iter > 0) {
    // starting EM algorithm to find values close to the MLE
    em_params = lognormal_mixture_em(em_iter, G, t, delta, X, better_initial_values, N_em, Niter_em, true, global_rng);
  } else if(show_output) {
    Rcout << "Skipping EM Algorithm" << "\n";
  }
  
  for (int iter = 0; iter < Niter; iter++) {
    // Starting empty objects for Gibbs Sampler
    if (iter == 0) {
      first_iter_gibbs(em_params, eta, beta, phi, em_iter, G, y, e0, sd, groups, X, use_W, delta, global_rng);
      
      // defining values for sintonizing the variance
      // of e0 proposal
      cte = 1.0;
      n_aceite = 0;
    }
    
    sd = 1.0 / sqrt(phi);
    
    // Data augmentation
    y_aug = augment(G, y, groups, delta, sd, beta, X, global_rng); 
    update_groups_gibbs(iter, use_W, em_params, G, y_aug, eta, beta, phi, X, groups, global_rng);
    
    // Computing number of observations allocated at each class
    n_groups = groups_table(G, groups);
    
    // ensuring that every class have, at least, 5 observations
    avoid_group_with_zero_allocation(n_groups, groups, G, N, global_rng);
    
    // updating all parameters
    update_gibbs_parameters(G, X, y_aug, n_groups, groups, eta, beta, phi, e0, global_rng);
    
    // this one is done separated from Gibbs because it's done via Metropolis-Hastings
    update_e0(G, iter, prop_aceite, n_aceite, cte, a, eta, e0, count, global_rng);
    
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
    
    newRow = arma::join_rows(beta.row(0),
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

/* SPARSE ALTERNATIVES */

// Function used to sample the latent groups for each observation.
arma::ivec sample_groups_sparse(const int& G, const arma::vec& y, const arma::vec& eta, 
                                const arma::vec& phi, const arma::mat& beta,
                                const arma::sp_mat& X, gsl_rng* rng_device,
                                const arma::ivec& groups_old) {
  int n = X.n_rows;
  int p = X.n_cols;
  arma::ivec vec_groups = groups_old;
  arma::ivec sequence = seq(0, G - 1);
  arma::uvec indexg;
  arma::mat Vg0 = arma::diagmat(repl(30.0, p));
  arma::mat Vg0_inv = arma::inv(Vg0);
  arma::vec xi;
  arma::rowvec xit;
  arma::sp_mat Xg;
  arma::sp_mat Xgt;
  arma::sp_mat Xgt_Xg;
  arma::vec yg;
  arma::mat S_inv(p, p);
  arma::mat S(p, p);
  arma::vec Xgt_yg;
  arma::vec Mg;
  arma::mat denom_mat(n, G);
  arma::vec probsEqual(G);
  probsEqual.fill(1.0 / G);
  
  double m;
  double yi;
  double a;
  double sigma2;
  double denom;
  
  for (int g = 0; g < G; g++) {
    indexg = arma::find(vec_groups == g);
    
    if(indexg.empty()) {
      continue;
    }
    
    Xg = arma::conv_to<arma::sp_mat>::from(arma::mat(X).rows(indexg));
    Xgt = Xg.t();
    yg = y(indexg);
    Xgt_yg = Xgt * yg;
    S = phi(g) * Xg.t() * Xg + Vg0_inv;
    if(arma::det(S) == 0) {
      S_inv = Vg0;
    } else {
      S_inv = arma::inv(phi(g) * Xg.t() * Xg + Vg0_inv,
                        arma::inv_opts::allow_approx);
    }
    
    Mg = X * S_inv * Xgt_yg;
    
    for (int i = 0; i < n; i++) {
      xit = X.row(i);
      yi = y(i);
      m = 1.0 / arma::as_scalar(xit * S_inv * xit.t());
      a = (m - phi(g)) * (vec_groups(i) == g) +
        m * (vec_groups(i) != g);
      
      sigma2 = (a + phi(g))/(a * phi(g));
      denom_mat(i, g) = 
        eta(g) * (R::dnorm(yi,
                  square(phi(g)) * sigma2 * (Mg(i) - yi / (a + phi(g))),
                  sqrt(sigma2), 
                  false) * (vec_groups(i) == g) +
                    R::dnorm(yi,
                             square(phi(g)) * sigma2 * ((a + phi(g)) * Mg(i)/(a + phi(g) + 1.0)),
                             sqrt(sigma2),
                             false) * (vec_groups(i) != g));
      
      if(g == (G - 1)) {
        denom = arma::sum(denom_mat.row(i));
        if(denom > 0) {
          vec_groups(i) = numeric_sample(sequence, denom_mat.row(i).t() / denom, 
                     rng_device);
        } else {
          vec_groups(i) = numeric_sample(sequence, probsEqual, rng_device);
        }
      }
    }
  }
  
  return(vec_groups);
}

// Function used to simulate survival time for censored observations.
// Here, delta is a vector such that delta(i) == 0 implies that the 
// i-th observation is censored. Otherwise, the i-th observation is a
// failure time.
arma::vec augment_sparse(const int& G, const arma::vec& y, const arma::ivec& groups,
                         const arma::ivec& delta, const arma::vec& sd, 
                         const arma::mat& beta, const arma::sp_mat& X, gsl_rng* rng_device) {
  
  int n = X.n_rows;
  arma::vec out = y;
  arma::mat mean = X * beta.t();
  
  arma::uvec censored_indexes = arma::find(delta == 0); // finding which observations are censored
  
  for (int i : censored_indexes) {
    out(i) = augment_yi(y(i), arma::as_scalar(mean(i, groups(i))), sd(groups(i)), rng_device);
  }
  
  return out;
}

// ** Auxiliary functions for EM algorithm **

// Compute weights matrix
arma::mat compute_W_sparse(const arma::vec& y, const arma::sp_mat& X, const arma::vec& eta, 
                           const arma::mat& beta, const arma::vec& sigma, 
                           const int& G) {
  
  int n = X.n_rows;
  double denom;
  arma::mat out(n, G);
  arma::mat mat_denom(n, G);
  arma::rowvec repl_vec = repl(1.0 / G, G).t();
  
  for(int g = 0; g < G; g++) {
    mat_denom.col(g) = eta(g) * arma::normpdf(y,
                  X * beta.row(g).t(),
                  repl(sigma(g), n));
  }
  
  for(int i = 0; i < n; i++) {
    denom = arma::sum(mat_denom.row(i));
    
    if(denom > 0) {
      out.row(i) = mat_denom.row(i) / denom;
    } else {
      out.row(i) = repl_vec;
    }
  }
  
  return(out);
}

// Create the latent variable z for censored observations
arma::vec augment_em_sparse(const arma::vec& y, const arma::ivec& delta,
                            const arma::sp_mat& X, const arma::mat& beta,
                            const arma::vec& sigma, const arma::mat& W,
                            const int& G) {
  int n = X.n_rows;
  arma::vec out = y;
  arma::mat mean = X * beta.t();
  arma::mat alpha_mat(n, G);
  arma::uvec censored_indexes = arma::find(delta == 0); // finding which observations are censored
  
  for(int g = 0; g < G; g++) {
    alpha_mat.col(g) = (y - mean.col(g))/sigma(g);
  }
  
  for (int i : censored_indexes) {
    out(i) = 0.0;
    
    for (int g = 0; g < G; g++) {
      out(i) += W(i, g) * compute_expected_value_truncnorm(arma::as_scalar(alpha_mat(i, g)), arma::as_scalar(mean(i, g)), sigma(g));;
    }
  }
  
  return out;
}

void update_beta_g_sparse(const arma::vec& colg, const arma::sp_mat& X, const int& g, const arma::vec& z, arma::mat& beta) {
  arma::sp_mat Wg;
  Wg = arma::diagmat(colg);
  
  if(arma::det(arma::mat(X.t() * Wg * X)) != 0.0) {
    beta.row(g) = arma::solve(arma::mat(X.t() * Wg * X),
             arma::vec(X.t() * Wg * z),
             arma::solve_opts::allow_ugly).t();
  }
}

void update_phi_g_sparse(const double& denom, const arma::uvec& censored_indexes, const arma::sp_mat& X, const arma::vec& colg, const arma::vec& y, const arma::vec& z,
                         const arma::vec& sd, const arma::mat& beta, const arma::vec& var, const int& g, const int& n, arma::vec& phi, gsl_rng* rng_device) {
  double alpha = 0.0;
  double quant = arma::as_scalar(arma::square(z - (X * beta.row(g).t())).t() * colg);
  
  for(int i : censored_indexes) {
    alpha = (y(i) - arma::as_scalar(X.row(i) * beta.row(g).t())) / sd(g);
    
    if(R::pnorm(alpha, 0.0, 1.0, true, false) < 1.0) {
      quant += colg(i) * var(g) * (1.0 - 
        (- alpha * R::dnorm(alpha, 0.0, 1.0, false))/(1.0 - R::pnorm(alpha, 0.0, 1.0, true, false)) -
        square(R::dnorm(alpha, 0.0, 1.0, false)/(1.0 - R::pnorm(alpha, 0.0, 1.0, true, false))));
    } else {
      quant += colg(i) * var(g) * (1.0 - 
        (-alpha * R::dnorm(alpha, 0.0, 1.0, false))/(1.0 - 0.999) -
        square(R::dnorm(alpha, 0.0, 1.0, false)/(1.0 - 0.999)));
    }
  }
  
  // to avoid numerical problems
  if (quant == 0) {
    phi(g) = rgamma_(0.5, 0.5, rng_device); // resample phi
  }
  
  // to avoid numerical problems
  if(phi(g) > 1e5 || phi.has_nan()) {
    phi(g) = rgamma_(0.5, 0.5, rng_device); // resample phi
  }
}

void update_em_parameters_sparse(const int& n, const int& G, arma::vec& eta, arma::mat& beta, arma::vec& phi, const arma::mat& W, const arma::sp_mat& X, const arma::vec& y,
                                 const arma::vec& z, const arma::ivec& delta, arma::vec& sd, gsl_rng* rng_device) {
  arma::vec colg(n);
  arma::vec var;
  double quant, denom, alpha;
  arma::uvec censored_indexes = arma::find(delta == 0);
  var = arma::square(sd);
  
  for (int g = 0; g < G; g++) {
    colg = W.col(g);
    
    eta(g) = arma::sum(colg) / n; // updating eta(g)
    update_beta_g_sparse(colg, X, g, z, beta); // updating beta for the group g
    update_phi_g_sparse(arma::sum(colg), censored_indexes, X, colg, y, z, sd, beta, var, g, n, phi, rng_device);
  }
}

// Compute model's log-likelihood to select the EM initial values
double loglik_em_sparse(const arma::vec& eta, const arma::mat& beta, const arma::vec& phi, const arma::mat& W, const arma::vec& z, const int& G,
                        const arma::sp_mat& X) {
  int N = X.n_rows;
  arma::vec sd = 1.0 / sqrt(phi);
  arma::mat mean = arma::mat(X * beta.t());
  double loglik = 0.0;
  
  for(int i = 0; i < N; i++) {
    for(int g = 0; g < G; g++) {
      loglik += W(i, g) * R::dnorm(z(i), arma::as_scalar(mean(i, g)),
                  sd(g), true);
    }
  }
  
  return loglik;
}

// EM for the lognormal mixture model.
arma::field<arma::mat> lognormal_mixture_em_sparse(const int& Niter, const int& G, const arma::vec& t, const arma::ivec& delta, const arma::sp_mat& X,
                                                   const bool& better_initial_values, const int& N_em,
                                                   const int& Niter_em, const bool& internal, gsl_rng* rng_device) {
  
  int n = X.n_rows;
  int k = X.n_cols;
  
  // initializing objects used on EM algorithm
  arma::vec y = log(t);
  arma::vec eta(G);
  arma::vec phi(G);
  arma::vec sd(G);
  arma::vec z(n);
  arma::mat W(n, G);
  arma::mat beta(G, k);
  arma::mat out(Niter, G * k + (G * 2));
  arma::field<arma::mat> out_internal_true(6);
  arma::field<arma::mat> out_internal_false(2);
  
  arma::field<arma::mat> em_params(6);
  arma::field<arma::mat> best_em(6);
  for(int iter = 0; iter < Niter; iter++) {
    if(iter == 0) { // sample starting values
      if(better_initial_values) {
        for (int init = 0; init < N_em; init ++) {
          em_params = lognormal_mixture_em_sparse(Niter_em, G, t, delta, X, false, 0, 0, true, rng_device);
          
          if(init == 0) {
            best_em = em_params;
            Rcout << "Initial LogLik: " << arma::as_scalar(best_em(5)) << "\n";
          } else {
            if(arma::as_scalar(em_params(5)) > arma::as_scalar(best_em(5))) { // comparing logliks
              Rcout << "Previous maximum: " << arma::as_scalar(best_em(5)) << " | New maximum: " << arma::as_scalar(em_params(5))  << "\n";
              best_em = em_params;
            }
          }
        }
        
        eta = best_em(0);
        beta = best_em(1);
        phi = best_em(2);
        W = best_em(3);
      } else {
        sample_initial_values_em(eta, phi, beta, sd, G, k, rng_device);
        W = compute_W_sparse(y, X, eta, beta, sd, G);
      }
    } else {
      sd = 1.0 / sqrt(phi);
      z = augment_em_sparse(y, delta, X, beta, sd, W, G);
      W = compute_W_sparse(z, X, eta, beta, sd, G);
      update_em_parameters_sparse(n, G, eta, beta, phi, W, X, y, z, delta, sd, rng_device);
    }
    
    // Fill the out matrix
    arma::rowvec newRow = 
      arma::join_rows(eta.row(0), 
                      beta.row(0),
                      phi.row(0));
    
    for (int g = 1; g < G; g++) {
      newRow = 
        arma::join_rows(newRow, 
                        eta.row(g), 
                        beta.row(g),
                        phi.row(g));
    }
    
    out.row(iter) = newRow;
  }
  
  if(internal) {
    out_internal_true(0) = eta;
    out_internal_true(1) = beta;
    out_internal_true(2) = phi;
    out_internal_true(3) = W;
    out_internal_true(4) = augment_em_sparse(y, delta, X, beta, 1.0 / sqrt(phi), W, G);
    out_internal_true(5) = loglik_em_sparse(eta, beta, phi, W, out_internal_true(4), G, X);
    
    return out_internal_true;
  } else {
    out_internal_false(0) = out;
    out_internal_false(1) = loglik_em_sparse(eta, beta, phi, W, augment_em_sparse(y, delta, X, beta, 1.0 / sqrt(phi), W, G), G, X);
    
    return out_internal_false;
  }
  
  return out_internal_false; // should never be reached
}

arma::mat lognormal_mixture_gibbs_implementation_sparse(const int& Niter, const int& em_iter,
                                                        const int& G, const arma::vec& t,
                                                        const arma::ivec& delta, 
                                                        const arma::sp_mat& X, const double& a, 
                                                        long long int starting_seed,
                                                        const bool& show_output, 
                                                        const int& chain_num,
                                                        const bool& use_W, const bool& better_initial_values,
                                                        const int& Niter_em, const int& N_em) {
  
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
  
  arma::vec y = log(t);
  
  // The output matrix should have Niter rows (1 row for each iteration) and
  // nColsOutput columns (1 column for each element).
  arma::mat out(Niter, nColsOutput);
  
  // The order of filling the output matrix matters a lot, since we can
  // make label switching accidentally. Latter this is going to be defined
  // so we can always fill the matrix in the correct order (by columns, always).
  
  arma::sp_mat Xt = X.t();
  arma::vec y_aug(N);
  arma::ivec n_groups(G);
  arma::mat means(N, G);
  arma::vec sd(G);
  arma::vec linearComb(N);
  arma::mat comb;
  arma::field<arma::mat> em_params(6);
  
  if(em_iter > 0) {
    // starting EM algorithm to find values close to the MLE
    em_params = lognormal_mixture_em_sparse(em_iter, G, t, delta, X, better_initial_values, N_em, Niter_em, true, global_rng);
  } else if(show_output) {
    Rcout << "Skipping EM Algorithm" << "\n";
  }
  
  // Starting other new values for MCMC algorithms
  arma::vec eta(G);
  arma::vec phi(G);
  arma::mat beta(G, p);
  arma::ivec groups(N);
  arma::ivec groups_start(N);
  
  arma::sp_mat Xg;
  arma::sp_mat Xgt;
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
  double u;
  arma::rowvec newRow;
  int m;
  int idx;
  
  for (int iter = 0; iter < Niter; iter++) {
    // Starting empty objects for Gibbs Sampler
    if (iter == 0) {
      if (em_iter != 0) {
        // we are going to start the values using the last EM iteration
        eta = em_params(0);
        beta = em_params(1);
        phi = em_params(2);
        
        groups_start = sample_groups_from_W(em_params(3));
      } else {
        eta = rdirichlet(repl(1, G), global_rng);
        
        for (int g = 0; g < G; g++) {
          phi(g) = rgamma_(0.5, 0.5, global_rng);
          beta.row(g) = rmvnorm(repl(0.0, p),
                   arma::diagmat(repl(50, p)),
                   global_rng).t();
        }
        
        sd = 1.0 / sqrt(phi);
        // Sampling classes for the observations
        groups_start = sample_groups_start(G, y, eta, global_rng);
      }
      
      // sampling value for e0
      e0 = rgamma_(1.0, 1.0, global_rng);
      
      // defining values for sintonizing the variance
      // of e0 proposal
      cte = 1;
      n_aceite = 0;
      
      if(use_W == false) {
        groups = sample_groups_sparse(G, y_aug, eta, phi, beta, X, global_rng, groups_start);
      } else {
        groups = groups_start;
      }
    }
    
    sd = 1.0 / sqrt(phi);
    
    // Data augmentation
    y_aug = augment_sparse(G, y, groups, delta, sd, beta, X, global_rng); 
    
    if (iter > 0) {
      if(use_W) {
        groups = sample_groups_from_W(em_params(3));
      } else {
        groups = sample_groups_sparse(G, y_aug, eta, phi, beta, X, global_rng, groups);
      }
    }
    
    // Computing number of observations allocated at each class
    n_groups = groups_table(G, groups);
    
    // ensuring that every class have, at least, 5 observations
    for(int g = 0; g < G; g++) {
      if(n_groups(g) == 0) {
        m = 0;
        while(m < 5) {
          idx = numeric_sample(seq(0, N),
                               repl(1.0 / N, N),
                               global_rng);
          
          if(n_groups(groups(idx)) > 5) {
            groups(idx) = g;
            m += 1;
          } 
        }
        
        // recalculating the number of groups
        n_groups = groups_table(G, groups);
      }
    }
    
    // Sampling new eta
    eta = rdirichlet(arma::conv_to<arma::Col<double>>::from(n_groups) + e0, 
                     global_rng);
    
    // For each g, sample new phi[g] and beta[g, _]
    for (int g = 0; g < G; g++) {
      indexg = arma::find(groups == g);
      Xg = arma::conv_to<arma::sp_mat>::from(arma::mat(X).rows(indexg));
      Xgt = Xg.t();
      yg = y_aug.rows(indexg);
      
      linearComb = yg - Xg * beta.row(g).t();
      
      // sampling new phi
      // the priori used was Gamma(0.01, 0.01)
      phi(g) = rgamma_(n_groups(g) / 2.0 + 0.01, 
          (1.0 / 2.0) * arma::as_scalar(linearComb.t() * linearComb) + 0.01,
          global_rng);
      
      // sampling beta new
      // the priori used was MNV(vec 0, diag 1000)
      comb = phi(g) * Xgt * Xg + arma::diagmat(repl(1.0 / 30.0, p));
      
      if(arma::det(comb) != 0) {
        Sg = arma::solve(comb,
                         arma::eye(X.n_cols, X.n_cols),
                         arma::solve_opts::allow_ugly);
        
        mg = arma::solve(comb,
                         phi(g) * Xgt * yg,
                         arma::solve_opts::allow_ugly);
        
        beta.row(g) = rmvnorm(mg, Sg, global_rng).t();
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
    
    b = cte * a;
    
    // updating the value of e0 (eta's dirichlet hyperparameter)
    e0_prop = rgamma_(b*e0, b, global_rng);
    
    log_eta_new = log(eta);
    
    log_alpha = arma::sum(log_eta_new) * (e0_prop - e0) + 9.0 * log(e0_prop/e0) -
      10.0*  G * (e0_prop - e0) + b * (e0_prop - e0) + (b * e0_prop - 1.0) *
      log(e0) - (b*e0 - 1.0) * log(e0_prop);
    
    u = runif_0_1(global_rng);
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
    
    newRow = arma::join_rows(beta.row(0),
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

struct GibbsWorker : public RcppParallel::Worker {
  const arma::vec& seeds; // starting seeds for each chain
  arma::cube& out; // store matrix iterations for each chain
  
  // other parameters used to fit the model
  const int& Niter;
  const int& em_iter;
  const int& G;
  const arma::vec& t;
  const arma::ivec& delta;
  const arma::mat& X;
  const double& a;
  const bool& show_output;
  const bool& use_W;
  const bool& better_initial_values;
  const int& N_em;
  const int& Niter_em;
  
  // Creating Worker
  GibbsWorker(const arma::vec& seeds, arma::cube& out, const int& Niter, const int& em_iter, const int& G, const arma::vec& t,
              const arma::ivec& delta, const arma::mat& X, const double& a, const bool& show_output, const bool& use_W, const bool& better_initial_values,
              const int& N_em, const int& Niter_em) :
    seeds(seeds), out(out), Niter(Niter), em_iter(em_iter), G(G), t(t), delta(delta), X(X), a(a), show_output(show_output), use_W(use_W), better_initial_values(better_initial_values), N_em(N_em), Niter_em(Niter_em) {}
  
  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t i = begin; i < end; ++i) {
      usleep(5000 * i); // avoid racing conditions
      out.slice(i) = lognormal_mixture_gibbs_implementation(Niter, em_iter, G, t, delta, X, a, seeds(i), show_output, i + 1, use_W, better_initial_values, Niter_em, N_em);
    }
  }
};

struct GibbsWorkerSparse : public RcppParallel::Worker {
  const arma::vec& seeds; // starting seeds for each chain
  arma::cube& out; // store matrix iterations for each chain
  
  // other parameters used to fit the model
  const int& Niter;
  const int& em_iter;
  const int& G;
  const arma::vec& t;
  const arma::ivec& delta;
  const arma::sp_mat& X;
  const double& a;
  const bool& show_output;
  const bool& use_W;
  const bool& better_initial_values;
  const int& N_em;
  const int& Niter_em;
  
  // Creating Worker
  GibbsWorkerSparse(const arma::vec& seeds, arma::cube& out, const int& Niter, const int& em_iter, const int& G, const arma::vec& t,
                    const arma::ivec& delta, const arma::sp_mat& X, const double& a, const bool& show_output, const bool& use_W, const bool& better_initial_values,
                    const int& N_em, const int& Niter_em) :
    seeds(seeds), out(out), Niter(Niter), em_iter(em_iter), G(G), t(t), delta(delta), X(X), a(a), show_output(show_output), use_W(use_W), better_initial_values(better_initial_values), N_em(N_em), Niter_em(Niter_em) {}
  
  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t i = begin; i < end; ++i) {
      usleep(5000 * i); // avoid racing conditions
      out.slice(i) = lognormal_mixture_gibbs_implementation_sparse(Niter, em_iter, G, t, delta, X, a, seeds(i), show_output, i + 1, use_W, better_initial_values, Niter_em, N_em);
    }
  }
};

// Function to call lognormal_mixture_gibbs_implementation with parallellization
// [[Rcpp::export]]
arma::cube lognormal_mixture_gibbs(const int& Niter, const int& em_iter, const int& G,
                                   const arma::vec& t, const arma::ivec& delta, 
                                   const arma::mat& X, const double& a, 
                                   const arma::vec& starting_seed, const bool& show_output,
                                   const int& n_chains, const bool& sparse, const bool& use_W,
                                   const bool& better_initial_values, const int& N_em, const int& Niter_em) {
  arma::cube out(Niter, (X.n_cols + 2) * G, n_chains); // initializing output object
  
  // Fitting in parallel
  if(sparse) {
    arma::sp_mat Y(X);
    GibbsWorkerSparse worker(starting_seed, out, Niter, em_iter, G, t, delta, Y, a, show_output, use_W, better_initial_values, N_em, Niter_em);
    RcppParallel::parallelFor(0, n_chains, worker);
  } else {
    GibbsWorker worker(starting_seed, out, Niter, em_iter, G, t, delta, X, a, show_output, use_W, better_initial_values, N_em, Niter_em);
    RcppParallel::parallelFor(0, n_chains, worker);
  }
  
  return out;
}

//[[Rcpp::export]]
arma::mat lognormal_mixture_em_implementation(const int& Niter, const int& G, const arma::vec& t,
                                              const arma::ivec& delta, const arma::mat& X, 
                                              long long int starting_seed, const bool& sparse,
                                              const bool& better_initial_values, const int& N_em,
                                              const int& Niter_em) {
  
  gsl_rng* global_rng = gsl_rng_alloc(gsl_rng_default);
  
  // setting global seed to start the sampler
  setSeed(starting_seed, global_rng);
  
  arma::field<arma::mat> out;
  if(sparse) {
    arma::sp_mat X_sp(X);
    out = lognormal_mixture_em_sparse(Niter, G, t, delta, X_sp, better_initial_values, N_em, Niter_em, false, global_rng);
  } else {
    out = lognormal_mixture_em(Niter, G, t, delta, X, better_initial_values, N_em, Niter_em, false, global_rng);
  }
  
  return out(0);
}
