# Function used to quick visualize the fitted values (survival estimate) on the data used to fit the model (via EM algorithm or Gibbs).

`plot_fit_on_data()` estimates survival/hazard for the data the model
was fitted on and plots the results.

## Usage

``` r
plot_fit_on_data(
  model,
  data,
  type = "survival",
  interval = "none",
  level = 0.95
)
```

## Arguments

- model:

  A `survival_ln_mixture` or `survival_ln_mixture_em` object.

- data:

  A [`data.frame()`](https://rdrr.io/r/base/data.frame.html) or
  `tibble()` containing the data used to fit the model. For appropriate
  behavior, should be the same object used to generate
  survival_ln_mixture/survival_ln_mixture_em objects.

- type:

  A character string specifying the type of plot. The default is
  "survival", but can be "hazard".

- interval:

  A character string specifying the type of interval to be plotted. The
  default is "none", but can be "credible". The EM algorithm does not
  provide confidence intervals and this parameter is only support for
  the Bayesian version (`survival_ln_mixture` object).

- level:

  A numeric value between 0 and 1 specifying the level of the confidence
  interval. The default is 0.95.

## Value

A list with two objects, one ggplot (`$ggplot`) with the predictions
plotted against the empirical data and a tibble with the predictions
(`$preds`).
