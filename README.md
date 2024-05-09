
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
initial values from a small run of the EM-Algorithm. The samplers are
implemented in `C++` using `RcppArmadillo` for the linear algebra
operations, `RcppGSL` for the random number generation and `OpenMP` for
parallellization.

## Dependencies

Before proceeding, make sure to have
[GSL](https://www.gnu.org/software/gsl/) installed on your system.
Below, there are specific tutorial on how to install these for each
operational system other than Windows (Windows users are probably fine
and ready to go).

The package also depends on OpenMP, which should be fine for both Linux
and Windows users. For MacOS users, there is a little guide bellow.

### Mac OS

Specifically, on Mac OS, running `brew install gsl libomp` on the
console/terminal should be enough for GSL and OpenMP. The OpenMP is
[known](https://github.com/Rdatatable/data.table) to have issues with
Mac’s default compiler (clang). Thus, there are some additional steps
that envolves modifying your **~/.R/Makevars**, appending

> CPPFLAGS += -I/opt/homebrew/opt/libomp/include -Xpreprocessor -fopenmp
>
> LDFLAGS += -L/opt/homebrew/opt/libomp/lib -lomp

to it.

If you don’t know how to do it or are having troubles accessing the
file, try opening the console/terminal and running `cat ~/.R/Makevars`.
If it’s empty or doesn’t exist, simply run

    mkdir ~/.R/ && touch ~/.R/Makevars && echo -e "LDFLAGS += -L/opt/homebrew/opt/libomp/lib -lomp\nCPPFLAGS += -I/opt/homebrew/opt/libomp/include -Xpreprocessor -fopenmp" >> ~/.R/Makevars

This will create the folder **~/.R**, create the file **Makevars**
inside this folder and append the flags cited above to it. Your Makevars
is something you don’t want to be messing around, so if after running
`cat ~/.R/Makevars` you see something, just run

    echo -e "LDFLAGS += -L/opt/homebrew/opt/libomp/lib -lomp\nCPPFLAGS += -I/opt/homebrew/opt/libomp/include -Xpreprocessor -fopenmp" >> ~/.R/Makevars

on your terminal to append the necessary flags to your **Makevars**.

After that, you’re ready to install the package.

## Installation

You can install the development version of `lnmixsurv` from
[GitHub](https://github.com/) with:

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
