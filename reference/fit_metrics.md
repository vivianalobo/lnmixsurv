# Function used to calculate some distance metrics between the predicted survival and the observed survival. `fit_metrics()` is used to calculate distance metrics between empirical and fitted survival for a predictions object, preferably the object `$preds` returned from `plot_fit_on_data()`.

Function used to calculate some distance metrics between the predicted
survival and the observed survival. `fit_metrics()` is used to calculate
distance metrics between empirical and fitted survival for a predictions
object, preferably the object `$preds` returned from
[`plot_fit_on_data()`](https://vivianalobo.github.io/lnmixsurv/reference/plot_fit_on_data.md).

## Usage

``` r
fit_metrics(preds, nobs = NULL, threshold = 0.005)
```

## Arguments

- preds:

  The `$preds` object from
  [`plot_fit_on_data()`](https://vivianalobo.github.io/lnmixsurv/reference/plot_fit_on_data.md)
  applied to the model. If not this one, should be a prediction `tibble`
  with the columns `time`, `strata` (if applicable), `estimate`,
  `.pred_survival`, `n.risk` (also `chain`, if necessary). It's
  important that the quantities `estimate` and `.pred_survival` are
  calculated for the same `time` and `strata`. It's highly recommended
  to simply use the object `$preds` returned from the function
  [`plot_fit_on_data()`](https://vivianalobo.github.io/lnmixsurv/reference/plot_fit_on_data.md).

- nobs:

  The number of observations used to fit the model. Can be ignored if
  `threshold` is set to 0. To easily calculate this value, use the
  function [`nobs()`](https://rdrr.io/r/stats/nobs.html) applied to the
  model object.

- threshold:

  Numeric value between 0 and 1. Times with `n.risk` below threshold \*
  nobs will be ignored. Default is 0.005 (0.5%). Important because the
  distance metrics may be too big if calculated in intervals without
  sufficient observations to be estimated.

## Value

A `tibble` with the following columns:

- `strata`: The stratas used to fit the model (if necessary).

- `n_strata`: The number of observations in the strata.

- `chain`: The chain of the Bayesian model (only if necessary).

- `metric`: Which metric is being calculated.

- `value`: Value for the metric.

For now, the following metrics are available and will be included:

- `MSE`: Mean Squared Error (the less the better).

- `MAE`: Mean Absolute Error (the less the better).

- `Hellinger Distance`: Hellinger distance, sometimes called Jeffreys
  distance (the less the better).

- `KS Distance`: Kolmogorov-Smirnov distance (the less the better).
