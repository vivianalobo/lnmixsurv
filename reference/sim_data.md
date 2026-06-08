# Simulated lognormal mixture data.

A simulated dataset with 10000 observations from a lognormal mixutre
model with 2 componentes.

## Usage

``` r
sim_data
```

## Format

### `sim_data`

A list with two componentes:

- \$data: A data frame with 10,000 rows and 3 columns:

  - y:

    observed survival time

  - delta:

    event indicator. 1 == event, 0 == censored.

  - x:

    binary covariate

- \$true_vals: A named vector with the true values used to generate the
  data.
