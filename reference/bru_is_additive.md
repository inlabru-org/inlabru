# Check for predictor expression additivity

Checks if a predictor expression is additive or not

## Usage

``` r
bru_is_additive(x, ...)

# S3 method for class 'character'
bru_is_additive(x, ..., verbose = FALSE)

# S3 method for class 'expression'
bru_is_additive(x, ..., verbose = FALSE)

# S3 method for class 'formula'
bru_is_additive(x, ..., verbose = FALSE)

# S3 method for class 'bru_pred_expr'
bru_is_additive(x, ...)

# S3 method for class 'bru_obs'
bru_is_additive(x, ...)

# S3 method for class 'bru_obs_list'
bru_is_additive(x, ...)
```

## Arguments

- x:

  A predictor `expression`, `formula`, or parse information
  `data.frame`.

- ...:

  Arguments passed on recursively.

- verbose:

  logical; if `TRUE`, print diagnostic parsing information.

## Value

`TRUE` if the expression is detected to be additive, `FALSE` otherwise.

## Examples

``` r
bru_is_additive(~ x + y)
#> [1] TRUE
bru_is_additive(~ x * y)
#> [1] FALSE
```
