# Predict from a Lognormal Mixture Model

Predict from a Lognormal Mixture Model

## Usage

``` r
# S3 method for class 'survival_ln_mixture'
predict(
  object,
  new_data,
  type,
  eval_time,
  interval = "none",
  level = 0.95,
  ...
)
```

## Arguments

- object:

  A `survival_ln_mixture` object.

- new_data:

  A data frame or matrix of new predictors.

- type:

  A single character. The type of predictions to generate. Valid options
  are:

  - `"time"` for the survival time. **not implmeented**

  - `"survival"` for the survival probability.

  - `"hazard"` for the hazard.

- eval_time:

  For type = "hazard" or type = "survival", the times for the
  distribution.

- interval:

  should interval estimates be added? Options are "none" and "credible".

- level:

  the tail area of the intervals. Default value is 0.95.

- ...:

  Not used, but required for extensibility.

## Value

A tibble of predictions. The number of rows in the tibble is guaranteed
to be the same as the number of rows in `new_data`.

## Note

Categorical predictors must be converted to factors before the fit,
otherwise the predictions will fail.

## Examples

``` r

# Categorical variables must be converted to factor before the fit.

require(survival)
# Wrong way of doing
set.seed(1)
mod <- survival_ln_mixture(Surv(time, status == 2) ~ factor(sex), lung, intercept = TRUE)

if (FALSE) { # \dontrun{
# this piece of code will throw error
predict(mod, data.frame(sex = 1), type = "survival", eval_time = 100)
} # }

# Correct way
lung$sex <- factor(lung$sex) # converting to factor before
set.seed(1)
mod2 <- survival_ln_mixture(Surv(time, status == 2) ~ sex, lung, intercept = TRUE)
# Note: the categorical predictors must be a character.
predict(mod2, data.frame(sex = "1"), type = "survival", eval_time = 100)
#> # A tibble: 1 × 2
#>   .pred            strata
#>   <list>           <fct> 
#> 1 <tibble [1 × 2]> sex=1 
```
