
<!-- README.md is generated from README.Rmd. Please edit that file -->

# persistencia

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

The goal of persistencia is to provide an easy interface to the bayesian
lognormal mixture model described in (incluir artigo da Viviana). An
usual formula-type model is implemented in `survival_ln_mixture`, with
the usual `suvival::Surv(time, event) ~ cov` interface. The model tries
to follow the [conventions for R modeling
packages](https://tidymodels.github.io/model-implementation-principles/),
and uses the [hardhat](https://hardhat.tidymodels.org/) structure.

An extension to the models defined by
[parsnip](https://parsnip.tidymodels.org/index.html) and
[censored](https://censored.tidymodels.org/articles/examples.html) is
also provided, adding the `survival_ln_mixture` engine to the
`survival_reg` model.

The underlying algorithm implementation is a gibbs sampler (incluir
detalhes sobre prioris e etc) and is implemmented in `C++`, using
`RcppArmadillo` for the linear algebra operations and `OpenMP` to
provide a parallelized way of generating multiple chains.

## Installation

You can install the development version of persistencia from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("vitorcapdeville/persistencia")
```

## Example

The model can be used with the usual formula interface or using the
`tidymodels` and `censored` structure.

``` r
library(persistencia)
library(dplyr)
lung$sex <- factor(lung$sex)
set.seed(1)
mod <- survival_ln_mixture(Surv(time, status == 2) ~ sex, lung, intercept = TRUE)
mod
#> survival_ln_mixture
#>  formula:      Surv(time, status == 2) ~ sex
#>  observations: 228
#>  predictors:   2
#> ------
#>             estimate  std.error
#> beta_a[1]  5.5938612 0.09095261
#> beta_a[2]  0.4000306 0.13214121
#> beta_b[1]  2.5266175 0.08858586
#> beta_b[2] -0.8485575 0.21621294
#> 
#> Auxiliary parameter(s):
#>           estimate   std.error
#> phi_a    1.4653951  0.23724308
#> phi_b   38.4461391 42.95998041
#> theta_a  0.9585784  0.01824815
```

For the tidymodels approach, when using `recipes`, some caution must be
taken due to the fact that recipes do not support inline functions, and
the original variables can only have their roles updated with
`update_role`.

``` r
library(tidymodels)
tidymodels_prefer()

dados <- lung |> mutate(outcome = Surv(time, status == 2))
rec <- recipe(outcome ~ sex, dados) |>
  step_mutate(sex = factor(sex)) |>
  step_dummy(sex)

mod_spec <- survival_reg() |>
  set_engine("survival_ln_mixture") |>
  set_mode("censored regression")

wflw <- workflow() |>
  add_model(mod_spec) |>
  add_recipe(rec)

set.seed(1)
mod_tidy <- wflw |> fit(dados)
mod_tidy
#> ══ Workflow [trained] ══════════════════════════════════════════════════════════
#> Preprocessor: Recipe
#> Model: survival_reg()
#> 
#> ── Preprocessor ────────────────────────────────────────────────────────────────
#> 2 Recipe Steps
#> 
#> • step_mutate()
#> • step_dummy()
#> 
#> ── Model ───────────────────────────────────────────────────────────────────────
#> survival_ln_mixture
#>  formula:      x ~ .
#>  observations: 228
#>  predictors:   2
#> ------
#>             estimate  std.error
#> beta_a[1]  5.5938612 0.09095261
#> beta_a[2]  0.4000306 0.13214121
#> beta_b[1]  2.5266175 0.08858586
#> beta_b[2] -0.8485575 0.21621294
#> 
#> Auxiliary parameter(s):
#>           estimate   std.error
#> phi_a    1.4653951  0.23724308
#> phi_b   38.4461391 42.95998041
#> theta_a  0.9585784  0.01824815
```

The predicitons are easy to obtain from a fit.

``` r
library(ggplot2)
library(purrr)
library(dplyr)
library(tidyr)

models <- list(formula = mod, tidymodels = mod_tidy)
new_data <- lung |> distinct(sex)
pred_sob <- map(models, ~ predict(.x, new_data, type = "survival", time = seq(120)))

bind_rows(pred_sob, .id = "modelo") %>%
  group_by(modelo) %>%
  dplyr::mutate(id = new_data$sex) %>%
  ungroup() %>%
  tidyr::unnest(cols = .pred) %>%
  ggplot(aes(x = .time, y = .pred_survival, col = id)) +
  geom_line() +
  theme_bw() +
  facet_wrap(~modelo)
```

<img src="man/figures/README-prediction-1.png" width="100%" />
