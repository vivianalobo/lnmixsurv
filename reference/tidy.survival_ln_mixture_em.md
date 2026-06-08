# Tidying method for a Lognormal Mixture model (fitted via Expectation-Maximization algorithm).

These method tidy the estimates from `survival_ln_mixture` fits into a
short summary. It doesn't contain uncertainty estimates since it's a
likelihood maximization algorithm.

## Usage

``` r
# S3 method for class 'survival_ln_mixture_em'
tidy(x, effects = "fixed", digits = NULL, ...)
```

## Arguments

- x:

  Fitted model object (survival_ln_mixture_em).

- effects:

  A character vector including one or more of `"fixed"` and
  `"auxiliary`.

- digits:

  How many significant digits should be displayed?

- ...:

  Not used.

## Value

A `data.frame` without rownames. When `effects="fixed"` (the default),
tidy.survival_ln_mixutre returns one row for each coefficient for each
component of the mixture with two columns:

- term:

  The name of the corresponding term in the model.

- estimate:

  A point estimate of the coefficient (last iteration value).

Setting `effects="auxiliary"` will select the precision and proportion
of mixture components parameters.

## Examples

``` r

require(survival)
lung$sex <- factor(lung$sex)
set.seed(1)
mod2 <- survival_ln_mixture_em(Surv(time, status == 2) ~ sex, lung)
tidy(mod2)
#> # A tibble: 4 × 2
#>   term          estimate
#>   <chr>            <dbl>
#> 1 (Intercept)_1    4.11 
#> 2 sex2_1           0.941
#> 3 (Intercept)_2    5.74 
#> 4 sex2_2           0.374
tidy(mod2, effects = c("fixed", "auxiliary"))
#> # A tibble: 8 × 2
#>   term          estimate
#>   <chr>            <dbl>
#> 1 (Intercept)_1    4.11 
#> 2 sex2_1           0.941
#> 3 (Intercept)_2    5.74 
#> 4 sex2_2           0.374
#> 5 phi_1            0.676
#> 6 phi_2            2.25 
#> 7 eta_1            0.199
#> 8 eta_2            0.801
```
