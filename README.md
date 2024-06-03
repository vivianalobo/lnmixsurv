
<!-- README.md is generated from README.Rmd. Please edit that file -->

# lnmixsurv

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![Codecov test
coverage](https://codecov.io/gh/vitorcapdeville/persistencia/branch/master/graph/badge.svg)](https://app.codecov.io/gh/vitorcapdeville/persistencia?branch=master)
[![R-CMD-check](https://github.com/vivianalobo/lnmixsurv/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/vivianalobo/lnmixsurv/actions/workflows/R-CMD-check.yaml)

<!-- badges: end -->

The `lnmixsurv` package provides an easy interface to the Bayesian
lognormal mixture model proposed by [Lobo, Fonseca and Alves,
2023](https://www.cambridge.org/core/journals/annals-of-actuarial-science/article/abs/lapse-risk-modeling-in-insurance-a-bayesian-mixture-approach/EDA511D313959D9A4040C51289A29B4A).

An usual formula-type model is implemented in `survival_ln_mixture`,
with the usual `suvival::Surv()` interface. The model tries to follow
the [conventions for R modeling
packages](https://tidymodels.github.io/model-implementation-principles/),
and uses the [hardhat](https://hardhat.tidymodels.org/) structure.

The underlying algorithm implementation is a Gibbs sampler which takes
initial values from a small run of the EM-Algorithm, with initial values
selection based on the log-likelihood. Besides the Bayesian approach,
the Expectation-Maximization approach (which focus on maximizing the
likelihood) for censored data is also available. The methods are
implemented in `C++` using `RcppArmadillo` for the linear algebra
operations, `RcppGSL` for the random number generation and seed control
and `RcppParallel` (since version 3.0.0) for parallelization.

## Dependencies

The only dependency is on GSL, so, make sure you have
[GSL](https://www.gnu.org/software/gsl/) installed before proceeding
Below, there are some basic guides on how to install these for each
operational system other than Windows (Windows users are probably fine
and ready to go).

### Mac OS

`brew install gsl` on the console/terminal should be enough for GSL.

### Linux

The installation of GSL on Linux is distro-specific. For the main
distros out-there:

- Debian/Ubuntu: `sudo apt-get install libgsl-dev`
- Arch: `sudo pacman -S gsl`
- Fedora: `sudo dnf install gsl-devel`

## Installation

You can install the latest development version of `lnmixsurv` from
[GitHub](https://github.com/):

``` r
# install.packages("devtools")
devtools::install_github("vivianalobo/lnmixsurv")
```

## parsnip and censored extension

An extension to the models defined by
[parsnip](https://parsnip.tidymodels.org/index.html) and
[censored](https://censored.tidymodels.org/articles/examples.html) is
also provided, adding the `survival_ln_mixture` engine to the
`parsnip::survival_reg()` model.

The following models, engines, and prediction type are
available/extended through `persistencia`:

| model        | engine                 | time | survival | linear_pred | raw | quantile | hazard |
|:-------------|:-----------------------|:-----|:---------|:------------|:----|:---------|:-------|
| survival_reg | survival_ln_mixture    | ✖    | ✔        | ✖           | ✖   | ✖        | ✔      |
| survival_reg | survival_ln_mixture_em | ✖    | ✔        | ✖           | ✖   | ✖        | ✔      |
