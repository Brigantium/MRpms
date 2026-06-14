# Print useful information of an `MRpms_conf` object

A generic function to summarize useful information of an `MRpms_conf`
object, which are computed from a `ConfPMS` call.

## Usage

``` r
# S3 method for class 'MRpms_conf'
summary(object, ...)
```

## Arguments

- object:

  An `MRpms_conf` object.

- ...:

  Other arguments. Currently, nothing more is avaliable.

## Value

The summarize of the `MRpms_conf` object, plus estimated uniform error
if computed and quantiles of pointwise errors if computed.
