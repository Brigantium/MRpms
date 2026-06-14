# Compute the volume of the uniform prediction set.

Compute the volume of the uniform prediction set for a fitted model.
First compute the modes in a given mesh using cross-validation, compute
the quantile for a given confidence level and later use Monte-Carlo
integration. Currently, only works in the univariate case.

## Usage

``` r
hP(
  X,
  Y,
  malla,
  dim = 1,
  h1,
  h2,
  p = floor(-log(eps, base = 10)),
  eps = 10^(-p),
  n = length(Y),
  conf.level = 0.95,
  k = attr(malla, "k"),
  lensopx = max(X) - min(X)
)
```

## Arguments

- X:

  Matrix with the coordinates of the covariable vector.

- Y:

  Vector with the response variable values.

- malla:

  A `malla` object from the `mallador` function.

- dim:

  Dimension of covariable vector.

- h1:

  Smooth parameter for covariable.

- h2:

  Smooth parameter for response variable.

- p:

  Number of decimal precision.

- eps:

  Tolerance for convergence in the partial Mean-Shift.

- n:

  Sample size.

- conf.level:

  Confidence level for the uniform prediction set.

- k:

  Number of initial guesses in every covariable point in the mesh.

- lensopx:

  Size of the support of covariable.

## Value

The volume of the estimated uniform prediction set.
