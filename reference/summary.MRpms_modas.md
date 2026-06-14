# Print useful information of an `MRpms_modas` object

A generic function to summarize useful information of an `MRpms_modas`
object, which are computed from a `PMS` call or, a more bare-bones
version, `PMSc` call.

## Usage

``` r
# S3 method for class 'MRpms_modas'
summary(object, ...)
```

## Arguments

- object:

  An `MRpms_modas` object.

- ...:

  Other arguments. Currently, nothing more is avaliable.

## Value

The maximum, minimum and mean number of modes computed, with certain
tolerance to avoid atipic values.
