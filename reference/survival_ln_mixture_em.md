# Lognormal mixture model - Expectation-Maximization Algorithm

`survival_ln_mixture_em()` fits an EM algorithm, as described in LOBO,
Viviana GR; FONSECA, Thaís CO; ALVES, Mariane B. Lapse risk modeling in
insurance: a Bayesian mixture approach. Annals of Actuarial Science, v.
18, n. 1, p. 126-151, 2024, for modelling mixtures of lognormal
distributions applied to survival data.

## Usage

``` r
survival_ln_mixture_em(
  formula,
  data,
  intercept = TRUE,
  iter = 50,
  mixture_components = 2,
  starting_seed = sample(1:2^28, 1),
  number_em_search = 200,
  iteration_em_search = 1,
  show_progress = FALSE,
  ...
)

# Default S3 method
survival_ln_mixture_em(formula, ...)

# S3 method for class 'formula'
survival_ln_mixture_em(formula, data, intercept = TRUE, ...)
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

  A positive integer specifying the number of iterations for the EM
  algorithm.

- mixture_components:

  number of mixture componentes \>= 2.

- starting_seed:

  Starting seed for the algorithm. If not specified by the user, uses a
  random integer between 1 and 2^28 This way we ensure, when the user
  sets a seed in R, that this is passed into the C++ code.

- number_em_search:

  Number of different EM's to search for maximum likelihoods.
  Recommended to leave, at least, at 100.

- iteration_em_search:

  Number of iterations for each of the EM's used to find the maximum
  likelihoods. Recommended to leave at small values, such as from 1 to
  5.

- show_progress:

  A logical. Should the progress of the EM algorithm be shown?

- ...:

  Not currently used, but required for extensibility.

## Value

An object of class `survival_ln_mixture_em` containing the following
elements:

- `em_iterations`: A data frame containing the EM iterations.

- `nobs`: The number of observations.

- `predictors_name`: The names of the predictors.

- `logLik`: The log-likelihood of the model.

- `mixture_groups`: The number of mixture groups.

- `blueprint`: The blueprint used to process the formula
