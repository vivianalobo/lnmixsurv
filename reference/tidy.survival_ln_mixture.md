# Tidying method for a Lognormal Mixture model.

These method tidy the estimates from `survival_ln_mixture` fits into a
summary.

## Usage

``` r
# S3 method for class 'survival_ln_mixture'
tidy(
  x,
  effects = "fixed",
  conf.int = FALSE,
  conf.level = 0.95,
  digits = NULL,
  ...
)
```

## Arguments

- x:

  Fitted model object (survival_ln_mixture).

- effects:

  A character vector including one or more of `"fixed"` and
  `"auxiliary`.

- conf.int:

  If `TRUE` columns for lower (`cred.low`) and upper (`cred.high`)
  bounds of the posterior uncertainty intervals are included.

- conf.level:

  A number between 0 and 1 indicating the desired probability mass to
  include in the intervals. Only used if `conf.int = TRUE`.

- digits:

  How many significant digits should be displayed?

- ...:

  Not used.

## Value

A `data.frame` without rownames. When `effects="fixed"` (the default),
tidy.survival_ln_mixutre returns one row for each coefficient for each
component of the mixture with three columns:

- term:

  The name of the corresponding term in the model.

- estimate:

  A point estimate of the coefficient (posterior median).

- std.error:

  A standard error for the point estimate based on
  [`mad`](https://rdrr.io/r/stats/mad.html). See the *Uncertainty
  estimates* section in `print.stanreg` for more details.

Setting `effects="auxiliary"` will select the precision and proportion
of mixture components parameters.

## Examples

``` r

require(survival)
lung$sex <- factor(lung$sex)
set.seed(1)
mod2 <- survival_ln_mixture(Surv(time, status == 2) ~ sex, lung)
tidy(mod2)
#> # A tibble: 4 × 3
#>   term          estimate std.error
#>   <chr>            <dbl>     <dbl>
#> 1 (Intercept)_1    5.76      0.108
#> 2 sex2_1           0.314     0.181
#> 3 (Intercept)_2    5.04      0.225
#> 4 sex2_2           0.822     0.345
tidy(mod2, conf.int = TRUE)
#> # A tibble: 4 × 5
#>   term          estimate std.error cred.low cred.high
#>   <chr>            <dbl>     <dbl>    <dbl>     <dbl>
#> 1 (Intercept)_1    5.76      0.108   5.52       5.99 
#> 2 sex2_1           0.314     0.181  -0.0698     0.690
#> 3 (Intercept)_2    5.04      0.225   4.59       5.48 
#> 4 sex2_2           0.822     0.345   0.179      1.55 
tidy(mod2, effects = c("fixed", "auxiliary"), conf.int = TRUE)
#> # A tibble: 7 × 5
#>   term          estimate std.error cred.low cred.high
#>   <chr>            <dbl>     <dbl>    <dbl>     <dbl>
#> 1 (Intercept)_1    5.76     0.108    5.52       5.99 
#> 2 sex2_1           0.314    0.181   -0.0698     0.690
#> 3 (Intercept)_2    5.04     0.225    4.59       5.48 
#> 4 sex2_2           0.822    0.345    0.179      1.55 
#> 5 phi_1            3.01     0.746    1.88       4.85 
#> 6 phi_2            0.550    0.115    0.363      0.801
#> 7 eta_1            0.512    0.0308   0.456      0.564
```
