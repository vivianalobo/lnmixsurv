# Augment data with information from a survival_ln_mixture_em object

Include information about hazard and survival distribution for each
individual in a dataset.

## Usage

``` r
# S3 method for class 'survival_ln_mixture_em'
augment(x, newdata, eval_time, ...)
```

## Arguments

- x:

  A `survival_ln_mixture_em` object.

- newdata:

  A [`base::data.frame()`](https://rdrr.io/r/base/data.frame.html) or
  `tibble::tiblle()` containing all the original predictors used to
  create x.

- eval_time:

  a vector with the times where the hazard and survival distribuition
  will be evaluated.

- ...:

  Not used.

## Value

A
[`tibble::tibble()`](https://tibble.tidyverse.org/reference/tibble.html)
with the original covariates and ther survvial and hazard distributions.
