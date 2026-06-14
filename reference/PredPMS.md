# Uniform prediction sets based on modal regression

Computes a uniform prediction set for the response value.

## Usage

``` r
PredPMS(
  data,
  modas,
  x.malla,
  h1 = 0.3,
  h2 = 0.5,
  dim.y = 1,
  eps = 1e-08,
  conf.level = 0.95,
  k = 10,
  len = 200,
  malla
)
```

## Arguments

- data:

  Matrix or data frame containing the sample. Currently, the last column
  is needed to be de response values.

- modas:

  List containing the estimated set of conditional local modes on each
  covariate point in `malla`. In case modas missing, PMSc is called to
  compute them. It must be a `MRpms_modas` class object.

- x.malla:

  Set of covariate points where the modes will be estimated. Not needed
  if `modas` is provided.

- h1:

  Smoothing parameter for covariates.

- h2:

  Smoothing parameter for response variable.

- dim.y:

  Response dimension.

- eps:

  Convergence tolerance and threshold to discriminate between modes.

- conf.level:

  Confidence level for the prediction set.

- k:

  Number of Y values per `x` point in the `malla` object created, if
  `malla` is not provided.

- len:

  Number of different `x` points in the `malla` object.

- malla:

  Matrix containing the initial points used to compute the modes using
  the Partial Mean-Shift algorithm. If not provided, `mallador` function
  is called. It must be a MRpms_malla class object.

## Value

An `MRpms_pred` object which contains,

- `modas`: an `MRpms_modas` object with the adjusted modes,

- `epsh`: the radius of the uniform prediction set.

- `conf.level`: the confidence level used to compute `epsh`.

Aditionally, the prediction set and the estimated modes are plotted in
case `dim.x` and `dim.y` are `1`.

## References

Chen, Y.-C., Genovese, C. R., Tibshirani, R. J. and Wasserman, L.
(2016). *Nonparametric modal regression*. The Annals of Statistics,
**44**(2), 489–514.

## Examples

``` r
pred <- PredPMS(twosines)

```
