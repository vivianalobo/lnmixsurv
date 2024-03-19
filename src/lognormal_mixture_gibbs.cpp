// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#include <omp.h>
#include <RcppArmadillo.h>
#include <RcppGSL.h>
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

double square(const double& x) {
  return x * x;
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
                         const arma::mat& X, gsl_rng* rng_device,
                         const arma::ivec& grupos_old) {
  int n = X.n_rows;
  int p = X.n_cols;
  arma::ivec vec_grupos = grupos_old;
  arma::ivec sequence = seq(0, G - 1);
  arma::uvec indexg;
  arma::mat Vg0 = arma::diagmat(repl(30.0, p));
  arma::mat Vg0_inv = arma::inv(Vg0);
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
    indexg = arma::find(vec_grupos == g);
    
    if(indexg.empty()) {
      continue;
    }
    
    Xg = X.rows(indexg);
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
      a = (m - phi(g)) * (vec_grupos(i) == g) +
        m * (vec_grupos(i) != g);
      
      sigma2 = (a + phi(g))/(a * phi(g));
      denom_mat(i, g) = 
        eta(g) * (R::dnorm(yi,
                  square(phi(g)) * sigma2 * (Mg(i) - yi / (a + phi(g))),
                  sqrt(sigma2), 
                  false) * (vec_grupos(i) == g) +
                    R::dnorm(yi,
                             square(phi(g)) * sigma2 * ((a + phi(g)) * Mg(i)/(a + phi(g) + 1.0)),
                             sqrt(sigma2),
                             false) * (vec_grupos(i) != g));
      
      if(g == (G - 1)) {
        denom = arma::sum(denom_mat.row(i));
        if(denom > 0) {
          vec_grupos(i) = numeric_sample(sequence, denom_mat.row(i).t() / denom, 
                     rng_device);
        } else {
          vec_grupos(i) = numeric_sample(sequence, probsEqual, rng_device);
        }
      }
    }
  }
  
  return(vec_grupos);
}

arma::ivec sample_groups_start(const int& G, const arma::vec& y, 
                               const arma::vec& eta, gsl_rng* rng_device) {
  int n = y.n_rows;
  arma::ivec vec_grupos(n);
  
  for (int i = 0; i < n; i++) {
    vec_grupos(i) = numeric_sample(seq(0, G - 1), eta, rng_device);
  }
  
  return(vec_grupos);
}

// Function used to simulate survival time for censored observations.
// Here, delta is a vector such that delta(i) == 0 implies that the 
// i-th observation is censored. Otherwise, the i-th observation is a
// failure time.

// This function can be vectorized with a little of work, which can
// improve performance. Take a look at
// https://cran.r-project.org/web/packages/RcppTN/RcppTN.pdf
// for sampling multiples numbers at the same time.
arma::vec augment(const int& G, const arma::vec& y, const arma::ivec& groups,
                  const arma::ivec& delta, const arma::vec& sd, 
                  const arma::mat& beta, const arma::mat& X, gsl_rng* rng_device) {
  
  int n = X.n_rows;
  arma::vec out(n);
  int g;
  double out_i;
  int count;
  arma::mat mean = X * beta.t();
  
  for (int i = 0; i < n; i++) {
    if (delta(i) == 1) {
      out(i) = y(i);
    } else {
      g = groups(i);
      out_i = y(i);
      count = 0;
      
      while(out_i <= y(i)) {
        out_i = rnorm_(arma::as_scalar(mean(i, g)),
                       sd(g), rng_device);
        
        // break if it seems like it's going to run forever
        // should never happen
        if(count > 10000) {
          out_i = y(i) + 0.01;
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
  
  return out;
}

/*
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

// ** Auxiliary functions for EM algorithm **
// These are NOT optimizated at all and we can make a lot of improvements here.
// Besides that, it converges usually with fewer than 200 iterations and runs
// significantly fast.

// Compute weights matrix
arma::mat compute_W(const arma::vec& y, const arma::vec& delta, 
                    const arma::mat& X, const arma::vec& eta, 
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
arma::vec augment_em(const arma::vec& y, const arma::vec& delta,
                     const arma::mat& X, const arma::mat& beta,
                     const arma::vec& sigma, const arma::mat& W,
                     const int& G) {
  int n = X.n_rows;
  arma::vec out = y;
  
  double quant;
  double alpha;
  double expected;
  arma::mat mean = X * beta.t();
  arma::mat alpha_mat(n, G);
  
  for(int g = 0; g < G; g++) {
    alpha_mat.col(g) = (y - mean.col(g))/sigma(g);
  }
  
  for (int i = 0; i < n; i++) {
    if(delta(i) == 0) {
      quant = 0;
      
      for (int g = 0; g < G; g++) {
        alpha = alpha_mat(i, g);
        
        if (R::pnorm(alpha, 0, 1, true, false) < 1) {
          expected = arma::as_scalar(mean(i, g)) + sigma(g) *
            (R::dnorm(alpha, 0, 1, false)/(1 - R::pnorm(alpha, 0, 1, true, false)));
        } else {
          expected = arma::as_scalar(mean(i, g)) + sigma(g) *
            (R::dnorm(alpha, 0, 1, false)/(1 - 0.999));
        }
        
        quant += W(i, g) * expected;
      }
      
      out(i) = quant;
    }
  }
  
  return(out);
}

arma::ivec sample_groups_from_W(const arma::mat& W) {
  int N = W.n_rows;
  arma::vec out(N);
  
  for(int r = 0; r < N; r++) {
    out(r) = W.row(r).index_max();
  }
  
  return(arma::conv_to<arma::ivec>::from(out));
}

// Internal implementation of the em_algorithm
arma::field<arma::mat> lognormal_mixture_em_internal(int Niter, int G,
                                                     arma::vec y, arma::vec delta,
                                                     arma::mat X, gsl_rng* rng_device) {
  arma::vec eta(G);
  int n = X.n_rows;
  int k = X.n_cols;
  
  arma::field<arma::mat> out(4);
  arma::mat beta(G, k);
  arma::vec phi(G);
  arma::vec sd(G);
  arma::vec var(G);
  arma::mat W(n, G);
  arma::sp_mat Wg(n, n);
  arma::vec z(n);
  arma::vec colg(n);
  
  double quant;
  double denom;
  double alpha;
  
  for(int iter = 0; iter < Niter; iter++) {
    if(iter == 0) { // sample starting values
      eta = rdirichlet(repl(1.0, G), rng_device);
      
      for (int g = 0; g < G; g++) {
        phi(g) = rgamma_(2.0, 8.0, rng_device);
        
        for (int c = 0; c < k; c++) {
          beta(g, c) = rnorm_(0.0, 15.0, rng_device);
        }
      }
      
      sd = 1.0 / sqrt(phi);
      W = compute_W(y, delta, X, eta, beta, sd, G);
    } else {
      sd = 1.0 / sqrt(phi);
      var = arma::square(sd);
      z = augment_em(y, delta, X, beta, sd, W, G);
      W = compute_W(z, delta, X, eta, beta, sd, G);
      
      for (int g = 0; g < G; g++) {
        colg = W.col(g);
        Wg = arma::diagmat(colg);
        
        eta(g) = arma::sum(colg) / n;
        
        beta.row(g) = arma::solve(X.t() * Wg * X,
                 X.t() * Wg * z,
                 arma::solve_opts::allow_ugly).t();
        
        quant = 0.0;
        denom = arma::sum(colg);
        
        for (int i = 0; i < n; i++) {
          quant += W(i, g) * square(z(i) - arma::as_scalar(X.row(i) * beta.row(g).t()));
          
          if (delta(i) == 0.0) {
            alpha = (y(i) - arma::as_scalar(X.row(i) * beta.row(g).t())) / sd(g);
            
            if(R::pnorm(alpha, 0.0, 1.0, true, false) < 1.0) {
              quant += W(i, g) * var(g) * (1.0 - 
                (- alpha * R::dnorm(alpha, 0, 1, false))/(1.0 - R::pnorm(alpha, 0, 1, true, false)) -
                square(R::dnorm(alpha, 0, 1, false)/(1.0 - R::pnorm(alpha, 0.0, 1.0, true, false))));
            } else {
              quant += W(i, g) * var(g) * (1.0 - 
                (-alpha * R::dnorm(alpha, 0, 1, false))/(1.0 - 0.999) -
                square(R::dnorm(alpha, 0.0, 1.0, false)/(1.0 - 0.999)));
            }
          }
        }
        
        phi(g) = denom / quant;
      }
    }
  }
  
  out(0) = eta;
  out(1) = beta;
  out(2) = phi;
  out(3) = W;
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
  
  arma::mat Xt = X.t();
  arma::vec y_aug(N);
  arma::ivec n_groups(G);
  arma::mat means(N, G);
  arma::vec sd(G);
  arma::vec linearComb(N);
  arma::mat comb;
  arma::field<arma::mat> em_params(4);
  
  if(em_iter > 0) {
    // starting EM algorithm to find values close to the MLE
    em_params = lognormal_mixture_em_internal(em_iter, G, y, arma::conv_to<arma::vec>::from(delta), X, global_rng);
  } else if(show_output) {
    Rcout << "Skipping EM Algorithm" << "\n";
  }
  
  // Starting other new values for MCMC algorithms
  arma::vec eta(G);
  arma::vec phi(G);
  arma::mat beta(G, p);
  arma::ivec groups(N);
  arma::ivec groups_start(N);
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
          phi(g) = rgamma_(2.0, 8.0, global_rng);
          beta.row(g) = rmvnorm(repl(0.0, p),
                   arma::diagmat(repl(7.0, p)),
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
      
      groups = sample_groups(G, y_aug, eta, phi, beta, X, global_rng, groups_start);
    }
    
    sd = 1.0 / sqrt(phi);
    
    // Data augmentation
    y_aug = augment(G, y, groups, delta, sd, beta, X, global_rng); 
    
    if (iter > 0) {
      groups = sample_groups(G, y_aug, eta, phi, beta, X, global_rng, groups);
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
      Xg = X.rows(indexg);
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

//[[Rcpp::export]]
arma::mat lognormal_mixture_em(int Niter, int G, arma::vec y, arma::vec delta,
                               arma::mat X, long long int starting_seed) {
  
  gsl_rng* global_rng = gsl_rng_alloc(gsl_rng_default);
  
  // setting global seed to start the sampler
  setSeed(starting_seed, global_rng);
  
  arma::vec eta(G);
  int n = X.n_rows;
  int k = X.n_cols;
  
  arma::mat out(Niter, G * k + (G * 2));
  arma::mat beta(G, k);
  arma::vec phi(G);
  arma::vec sd(G);
  arma::vec var(G);
  arma::mat W(n, G);
  arma::sp_mat Wg(n, n);
  arma::vec z(n);
  arma::vec colg(n);
  
  double quant;
  double denom;
  double alpha;
  
  for(int iter = 0; iter < Niter; iter++) {
    if(iter == 0) { // sample starting values
      eta = rdirichlet(repl(1.0, G), global_rng);
      
      for (int g = 0; g < G; g++) {
        phi(g) = rgamma_(2.0, 8.0, global_rng);
        
        for (int c = 0; c < k; c++) {
          beta(g, c) = rnorm_(0.0, 15.0, global_rng);
        }
      }
      
      sd = 1.0 / sqrt(phi);
      W = compute_W(y, delta, X, eta, beta, sd, G);
    } else {
      sd = 1.0 / sqrt(phi);
      var = arma::square(sd);
      z = augment_em(y, delta, X, beta, sd, W, G);
      W = compute_W(z, delta, X, eta, beta, sd, G);
      
      for (int g = 0; g < G; g++) {
        colg = W.col(g);
        Wg = arma::diagmat(colg);
        
        eta(g) = arma::sum(colg) / n;
        beta.row(g) = arma::solve(X.t() * Wg * X,
                 X.t() * Wg * z,
                 arma::solve_opts::allow_ugly).t();
        
        quant = 0.0;
        denom = arma::sum(colg);
        
        for (int i = 0; i < n; i++) {
          quant += W(i, g) * square(z(i) - arma::as_scalar(X.row(i) * beta.row(g).t()));
          
          if (delta(i) == 0.0) {
            alpha = (y(i) - arma::as_scalar(X.row(i) * beta.row(g).t())) / sd(g);
            
            if(R::pnorm(alpha, 0.0, 1.0, true, false) < 1.0) {
              quant += W(i, g) * var(g) * (1.0 - 
                (- alpha * R::dnorm(alpha, 0, 1, false))/(1.0 - R::pnorm(alpha, 0, 1, true, false)) -
                square(R::dnorm(alpha, 0, 1, false)/(1.0 - R::pnorm(alpha, 0.0, 1.0, true, false))));
            } else {
              quant += W(i, g) * var(g) * (1.0 - 
                (-alpha * R::dnorm(alpha, 0, 1, false))/(1.0 - 0.999) -
                square(R::dnorm(alpha, 0.0, 1.0, false)/(1.0 - 0.999)));
            }
          }
        }
        
        phi(g) = denom / quant;
      }
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
  return(out);
}

//[[Rcpp::export]]
arma::mat lognormal_mixture_sem(int Niter, int G, arma::vec y, arma::vec delta,
                                arma::mat X, long long int starting_seed) {
  
  gsl_rng* global_rng = gsl_rng_alloc(gsl_rng_default);
  
  // setting global seed to start the sampler
  setSeed(starting_seed, global_rng);
  
  arma::vec eta(G);
  int n = X.n_rows;
  int k = X.n_cols;
  
  arma::mat out(Niter, G * k + (G * 2));
  arma::mat beta(G, k);
  arma::vec phi(G);
  arma::vec sd(G);
  arma::vec var(G);
  arma::mat W(n, G);
  arma::ivec grupos(n);
  arma::mat W_grupos(n, G, arma::fill::zeros);
  arma::ivec n_grupos(G);
  arma::sp_mat Wg(n, n);
  arma::vec colg(n);
  arma::vec z(n);
  int idx;
  int m;
  
  double quant;
  double denom;
  double alpha;
  
  for(int iter = 0; iter < Niter; iter++) {
    if(iter == 0) { // sample starting values
      eta = rdirichlet(repl(1.0, G), global_rng);
      
      for (int g = 0; g < G; g++) {
        phi(g) = rgamma_(2.0, 8.0, global_rng);
        
        for (int c = 0; c < k; c++) {
          beta(g, c) = rnorm_(0.0, 15.0, global_rng);
        }
      }
      
      sd = 1.0 / sqrt(phi);
      W = compute_W(y, delta, X, eta, beta, sd, G);
      
      for(int i = 0; i < n; i++) {
        grupos(i) = numeric_sample(seq(0, G - 1), W.row(i).t(), global_rng);
      }
      
      // Computing number of observations allocated at each class
      n_grupos = groups_table(G, grupos);
      
      // ensuring that every class have, at least, 5 observations
      for(int g = 0; g < G; g++) {
        if(n_grupos(g) == 0) {
          m = 0;
          while(m < 5) {
            idx = numeric_sample(seq(0, n),
                                 repl(1.0 / n, n),
                                 global_rng);
            
            if(n_grupos(grupos(idx)) > 5) {
              grupos(idx) = g;
              m += 1;
            } 
          }
          
          // recalculating the number of groups
          n_grupos = groups_table(G, grupos);
        }
      }
      
      for(int i = 0; i < n; i++) {
        W_grupos(i, grupos(i)) = 1;
      }
    } else {
      sd = 1.0 / sqrt(phi);
      var = arma::square(sd);
      z = augment_em(y, delta, X, beta, sd, W_grupos, G);
      W = compute_W(z, delta, X, eta, beta, sd, G);
      
      W_grupos.fill(0);
      
      for(int i = 0; i < n; i++) {
        grupos(i) = numeric_sample(seq(0, G - 1), W.row(i).t(), global_rng);
      }
      
      // Computing number of observations allocated at each class
      n_grupos = groups_table(G, grupos);
      
      // ensuring that every class have, at least, 5 observations
      for(int g = 0; g < G; g++) {
        if(n_grupos(g) == 0) {
          m = 0;
          while(m < 5) {
            idx = numeric_sample(seq(0, n),
                                 repl(1.0 / n, n),
                                 global_rng);
            
            if(n_grupos(grupos(idx)) > 5) {
              grupos(idx) = g;
              m += 1;
            } 
          }
          
          // recalculating the number of groups
          n_grupos = groups_table(G, grupos);
        }
      }
      
      for(int i = 0; i < n; i++) {
        W_grupos(i, grupos(i)) = 1;
      }
      
      for (int g = 0; g < G; g++) {
        colg = W_grupos.col(g);
        Wg = arma::diagmat(colg);
        
        eta(g) = arma::sum(colg) / n;
        beta.row(g) = arma::solve(X.t() * Wg * X,
                 X.t() * Wg * z,
                 arma::solve_opts::allow_ugly).t();
        
        quant = 0.0;
        denom = arma::sum(colg);
        
        for (int i = 0; i < n; i++) {
          quant += W_grupos(i, g) * square(z(i) - arma::as_scalar(X.row(i) * beta.row(g).t()));
          
          if (delta(i) == 0.0) {
            alpha = (y(i) - arma::as_scalar(X.row(i) * beta.row(g).t()))/sd(g);
            
            if(R::pnorm(alpha, 0.0, 1.0, true, false) < 1.0) {
              quant += W_grupos(i, g) * var(g) * (1.0 - 
                (- alpha * R::dnorm(alpha, 0, 1, false))/(1.0 - R::pnorm(alpha, 0, 1, true, false)) -
                square(R::dnorm(alpha, 0, 1, false)/(1.0 - R::pnorm(alpha, 0.0, 1.0, true, false))));
            } else {
              quant += W_grupos(i, g) * var(g) * (1.0 - 
                (-alpha * R::dnorm(alpha, 0, 1, false))/(1.0 - 0.999) -
                square(R::dnorm(alpha, 0.0, 1.0, false)/(1.0 - 0.999)));
            }
          }
        }
        
        phi(g) = denom / quant;
      }
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
  
  return(out);
}