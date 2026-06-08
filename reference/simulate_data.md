# Function to simulate survival data from a mixture of normal distribution.

`simulate_data()` simulates data from a mixture model.

## Usage

``` r
simulate_data(
  n = 4000,
  mixture_components = 2,
  k = 2,
  percentage_censored = 0.4,
  starting_seed = sample(1:2^28, 1)
)
```

## Arguments

- n:

  Number of observations desired.

- mixture_components:

  Number of mixtures to include in the generation of the data.

- k:

  number of covariates generated (the total of covariates will be
  intercept + (k - 1) covariates).

- percentage_censored:

  Percentage of censored observations (defined as decimal value between
  0 and 1). This will generate a delta vector in which 1 is an event
  that ocurred and 0 is a censored observation.

- starting_seed:

  Seed to start the random number generation.

## Value

A list with two elements: `data` and `real_values`. The `data` element
is a tibble with the simulated data. The `real_values` is a tibble with
the real values of the parameters used to generate the data.
