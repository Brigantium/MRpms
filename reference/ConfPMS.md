# Pointwise and uniform confidence sets

Computes pointwise and uniform confidence sets for conditional local
modes and plots them along with the estimated modes over the given
sample. (Warning: a large number of bootstrap replications increases
execution time).

## Usage

``` r
ConfPMS(
  data,
  modas,
  malla,
  h1 = 0.3,
  h2 = 0.5,
  eps = 1e-08,
  dim.y = 1,
  k = 10,
  len = 200,
  conf.level = 0.95,
  B = 30,
  type = 2,
  seed = 2026
)
```

## Arguments

- data:

  Matrix or data frame containing the sample. Currently, the last column
  is needed to be the response values.

- modas:

  List containing the estimated set of conditional local modes on each
  covariate point in `malla`. In case `modas` missing, `PMSc`is called
  to compute them. It must be a MRpms_modas class object.

- malla:

  Matrix containing the initial points used to compute the modes using
  the Partial Mean-Shift algorithm. If not provided, `mallador`function
  is called. It must be a MRpms_malla class object.

- h1:

  Smoothing parameter for covariates.

- h2:

  Smoothing parameter for response variable.

- eps:

  PMS convergence tolerance.

- dim.y:

  Response dimension.

- k:

  Number of Y values per x point in the `malla` object created.

- len:

  Number of different X points in the `malla` object.

- conf.level:

  Confidence level for the sets.

- B:

  Number of *bootstrap* replications computed to estimate the confidence
  sets.

- type:

  Numerical argument. "0" means computing pointwise sets only, "1" for
  uniform sets only and "2" for both. Set to "2" by default.

- seed:

  Seed for the *bootstrap* simulations.

## Value

A list. First entry is an `MRpms_modas` object containing the estimated
set of conditional local modes for each `x` point in the `malla`. If
`type` is 0 or 2, returns a vector with the estimated pointwise errors
(`deltax`). If `type` is 1 or 2, returns the estimated uniform error
(`delta`). In case `dim.y` and `dim.x` are `1`, a plot with the desired
regions is printed.

## References

Chen, Y.-C., Genovese, C. R., Tibshirani, R. J. and Wasserman, L.
(2016). *Nonparametric modal regression*. The Annals of Statistics,
**44**(2), 489–514.

## Examples

``` r
conf <- ConfPMS(twosines, B = 10)


```
