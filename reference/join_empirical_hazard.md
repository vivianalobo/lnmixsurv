# Function used to join the empirical hazard to the data

`join_empirical_hazard()` takes a Kaplan Meier empirical estimate which
includes the survival estimates and joins the hazard estimates to it.

## Usage

``` r
join_empirical_hazard(km)
```

## Arguments

- km:

  Kaplan-Meier estimates, i.e., object generated after running
  broom::tidy(survfit_obj), in which survfit_obj is a survfit object.
  Can also be a survfit object.

## Value

The same object as inputed, but with the hazard estimates column
(hazard_estimate) joined to it.
