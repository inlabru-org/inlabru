# Extract a bru model fit into a tidy tibble

This function extracts the fixed effect coefficients or hyperparameters
from a fitted `bru` object and returns them in a tidy tibble format. See
[`generics::tidy()`](https://generics.r-lib.org/reference/tidy.html) for
more details on the tidy data format.

## Usage

``` r
# S3 method for class 'bru'
tidy(x, effects = "fixed", ...)
```

## Arguments

- x:

  A fitted `bru` object.

- effects:

  `"fixed"` (default) or `"hyperpar"`.

- ...:

  Unused.

## Value

A tibble with one row per term.
