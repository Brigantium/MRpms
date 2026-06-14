# Plot method for `MRpms_pred` class with dimension `1`.

Plot method for `MRpms_pred` objects, which are computed from a
`PredPMS` call and contains information to construct prediction
intervals for the computed modes.

## Usage

``` r
# S3 method for class 'MRpms_pred'
plot(x, data, ...)
```

## Arguments

- x:

  An `MRpms_pred` object.

- data:

  The data from which the modes were computed.

- ...:

  Other parameters

## Value

A plot with the `data`, in blank circles and the uniform prediction
intervals in gray.

## See also

[`PredPMS()`](https://brigantium.github.io/MRpms/reference/PredPMS.md)
