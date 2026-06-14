# Compute the cross-validation like function suggested by Zhou and Huang

This function mask a C function which preforms the cross validation
suggested by Zhou and Huang to search the optimal smooth parameter.

## Usage

``` r
CVC(
  data = cbind(X, Y),
  X = data[, 1:dim],
  Y = data[, dim + 1],
  malla = mallador(data),
  dim,
  h1,
  h2,
  p = floor(-log(eps, base = 10)),
  eps = 10^(-p),
  n = length(Y),
  k = attr(malla, "k")
)
```

## Arguments

- data:

  Matrix or data frame containing the sample.

- X:

  Matrix with the covariable points.

- Y:

  Vector with the response values.

- malla:

  An `MRpms_malla` object returned by the `mallador` function with the
  initial guesses to comppute the modes.

- dim:

  Dimension of covariable.

- h1:

  Smooth parameter for covariable.

- h2:

  Smooth parameter for response.

- p:

  Decimal precision.

- eps:

  Convergence tolerance.

- n:

  Sample size.

- k:

  Number of guesses for each covariable point in mesh.

## Value

The value of the cross-validation like function suggested by Zhou and
Huang,

\$\$\displaystyle \frac{1}{n}\sum\_{i=1}^n
d^2\left(\widehat{M}\_{n,\mathbf{h},-i}\left(X_i\right),Y_i\right)\widehat{N}\_{\mathbf{h},-i}^2\left(X_i\right)w\left(X_i\right),\$\$

where \\\widehat M\_{n,\mathbf{h},-i}(X_i)\\ is the estimated modes in
\\X_i\\ without consider \\X_i\\ and
\\\widehat{N}\_{\mathbf{h},-i}^2\left(X_i\right)w\left(X_i\right)\\ is
the number of modes estimated in \\X_i\\.

## References

Zhou, H. and Huang, X. (2019). *Bandwidth selection for nonparametric
modal regression*. Communication in Statistics - Simulation and
Computation, **48**(4), 968–984.
