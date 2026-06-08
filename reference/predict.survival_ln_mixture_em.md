# Predict from a lognormal_em Mixture Model fitted using EM algorithm.

Predict from a lognormal_em Mixture Model fitted using EM algorithm.

## Usage

``` r
# S3 method for class 'survival_ln_mixture_em'
predict(object, new_data, type, eval_time, ...)
```

## Arguments

- object:

  A `survival_ln_mixture_em` object.

- new_data:

  A data frame or matrix of new predictors.

- type:

  A single character. The type of predictions to generate. Valid options
  are:

  - `"survival"` for the survival probability.

  - `"hazard"` for the hazard theoretical hazard.

- eval_time:

  For type = "hazard" or type = "survival", the times for the
  distribution.

- ...:

  Not used, but required for extensibility.

## Value

A tibble of predictions. The number of rows in the tibble is guaranteed
to be the same as the number of rows in `new_data`.

## Note

Categorical predictors must be converted to factors before the fit,
otherwise the predictions will fail.
