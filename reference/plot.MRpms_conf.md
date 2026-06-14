# Plot method for `MRpms_conf` class with dimension `1`.

Plot method for `MRpms_conf` objects, which are computed from a
`ConfPMS` call and contains information to construct confidence
intervals for the computed modes.

## Usage

``` r
# S3 method for class 'MRpms_conf'
plot(x, data, ...)
```

## Arguments

- x:

  An `MRpms_conf` object.

- data:

  The data from which the modes were computed.

- ...:

  Other parameters

## Value

A plot with the `data`, in blank circles, the modes in red, the
pointwise confidence intervals in blue and the uniform confidence
intervals in gray.

## See also

[`ConfPMS()`](https://brigantium.github.io/MRpms/reference/ConfPMS.md)
