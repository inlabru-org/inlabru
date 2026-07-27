# Glance at a bru model fit

This function returns a one-row tibble of model summaries from a fitted
`bru` object. See
[`generics::glance()`](https://generics.r-lib.org/reference/glance.html)
for more details on the glance format.

## Usage

``` r
# S3 method for class 'bru'
glance(x, ...)
```

## Arguments

- x:

  A fitted `bru` object.

- ...:

  Unused.

## Value

A one-row tibble of model-level summaries. `elapsed` is wall-clock
seconds summed across the phases in `fit$bru_timings`, not CPU time.
`nobs` is the total row count summed across all likelihoods, returned as
`NA` for models that involve `family = "cp"`, as "observation count" is
not meaningful for such models.
