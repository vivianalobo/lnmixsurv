# Lognormal mixture model - Gibbs sampler

`survival_ln_mixture()` fits a Bayesian lognormal mixture model with
Gibbs sampling (optional EM algorithm to find local maximum at the
likelihood function), as described in LOBO, Viviana GR; FONSECA, Thaís
CO; ALVES, Mariane B. Lapse risk modeling in insurance: a Bayesian
mixture approach. Annals of Actuarial Science, v. 18, n. 1, p. 126-151,
2024.

## Usage

``` r
survival_ln_mixture(
  formula,
  data,
  intercept = TRUE,
  iter = 1000,
  warmup = floor(iter/10),
  thin = 1,
  chains = 1,
  cores = 1,
  mixture_components = 2,
  show_progress = FALSE,
  em_iter = 0,
  starting_seed = sample(1:2^28, 1),
  use_W = FALSE,
  number_em_search = 200,
  iteration_em_search = 1,
  fast_groups = TRUE,
  ...
)

# Default S3 method
survival_ln_mixture(formula, ...)

# S3 method for class 'formula'
survival_ln_mixture(formula, data, intercept = TRUE, ...)
```

## Arguments

- formula:

  A formula specifying the outcome terms on the left-hand side, and the
  predictor terms on the right-hand side. The outcome must be a
  [survival::Surv](https://rdrr.io/pkg/survival/man/Surv.html) object.

- data:

  A **data frame** containing both the predictors and the outcome.

- intercept:

  A logical. Should an intercept be included in the processed data?

- iter:

  A positive integer specifying the number of iterations for each chain
  (including warmup).

- warmup:

  A positive integer specifying the number of warmup (aka burnin)
  iterations per chain. The number of warmup iterations should be
  smaller than iter.

- thin:

  A positive integer specifying the period for saving samples.

- chains:

  A positive integer specifying the number of Markov chains.

- cores:

  A positive integer specifying the maximum number of cores to run the
  chains. Setting this to a value bigger than 1 will automatically
  trigger the parallel mode

- mixture_components:

  number of mixture componentes \>= 2.

- show_progress:

  Indicates if the code shows the progress of the EM algorithm and the
  Gibbs Sampler.

- em_iter:

  A positive integer specifying the number of iterations for the EM
  algorithm. The EM algorithm is performed before the Gibbs sampler to
  find better initial values for the chains. On simulations, values
  lower than 200 seems to work nice.

- starting_seed:

  Starting seed for the sampler. If not specified by the user, uses a
  random integer between 1 and 2^28 This way we ensure, when the user
  sets a seed in R, that this is passed into the C++ code.

- use_W:

  Specifies is the W (groups weight's matrix for each observation)
  should be used from EM. It holds W constant through the code,
  resulting in a faster Bayesian Inference (close to what Empirical
  Bayes would do). It may helps generating credible intervals for the
  survival and hazard curves, using the information from the previous EM
  iteration. Make sure the EM have converged before setting this
  parameter to true. In doubt, leave this as FALSE, the default.

- number_em_search:

  Number of different EM's to search for maximum likelihoods.
  Recommended to leave, at least, at 100. This value can be set to 0 to
  disable the search for maximum likelihood initial values.

- iteration_em_search:

  Number of iterations for each of the EM's used to find the maximum
  likelihoods. Recommended to leave at small values, such as from 1 to
  5.

- fast_groups:

  Use fast computation of groups allocations probabilities, defaults to
  TRUE. Setting it to FALSE can increase the computation time (a lot)
  but it's worth trying if the chains are not converging.

- ...:

  Not currently used, but required for extensibility.

## Value

A `survival_ln_mixture` object, which is a list with the following
componentes:

- posterior:

  A
  [posterior::draws_matrix](https://mc-stan.org/posterior/reference/draws_matrix.html)
  with the posterior of the parameters of the model.

- nobs:

  A integer holding the number of observations used to generate the fit.

- blueprint:

  The blueprint component of the output of
  [hardhat::mold](https://hardhat.tidymodels.org/reference/mold.html)

## Note

Categorical predictors must be converted to factors before the fit,
otherwise the predictions will fail.

## Examples

``` r

# Formula interface
library(survival)
set.seed(1)
mod <- survival_ln_mixture(Surv(time, status == 2) ~ NULL, lung, intercept = TRUE)
```
