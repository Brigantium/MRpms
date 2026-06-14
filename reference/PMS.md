# Preforms the Partial Mean-Shift Algorithm

Finds the local modes of `data` conditional on a given set of covariate
values.

## Usage

``` r
PMS(
  data,
  x.malla = NULL,
  h1 = 0.3,
  h2 = 0.5,
  eps = 1e-08,
  k = 10,
  len = 200,
  dim.y = 1,
  malla = NULL
)
```

## Arguments

- data:

  Matrix or data frame containing the sample. Currently, the last column
  is needed to be the response values.

- x.malla:

  Set of covariate points where the modes will be estimated.

- h1:

  Smoothing parameter for covariates.

- h2:

  Smoothing parameter for response variable.

- eps:

  Convergence tolerance.

- k:

  Number of Y values per x point in the `malla` object created.

- len:

  Number of different X points in the `malla` object.

- dim.y:

  Response dimension.

- malla:

  Matrix containing the initial points used to compute the modes with
  the Partial Mean-Shift algorithm. If not provided, `mallador`function
  is called. It must be a `MRpms_malla` class object.

## Value

A `MRpms_modas` object, which contains

- `$modas`: contains the estimated local modes conditional on the values
  of the covariate appearing in `x.malla`;

- `$x.malla`: with each covariable point where modes were computed;

- `$dims`: Covariable and response dimensions.

## References

Chen, Y.-C., Genovese, C. R., Tibshirani, R. J. and Wasserman, L.
(2016). *Nonparametric modal regression*. The Annals of Statistics,
**44**(2), 489–514.

## Examples

``` r
modas <- PMS(twosines)
plot(modas,twosines,pch = 19, col = "red")

```
