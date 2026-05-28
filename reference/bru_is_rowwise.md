# Check for predictor rowwise evaluability

Checks if a predictor expression may be evaluated rowwise

## Usage

``` r
bru_is_rowwise(x, ...)

# S3 method for class 'bru_pred_expr'
bru_is_rowwise(x, ...)

# S3 method for class 'bru_obs'
bru_is_rowwise(x, ...)

# S3 method for class 'bru_obs_list'
bru_is_rowwise(x, ...)

# S3 method for class 'bru_comp_list'
bru_is_rowwise(x, ...)

# S3 method for class 'bru_comp'
bru_is_rowwise(x, ...)

# S3 method for class 'bru_mapper'
bru_is_rowwise(x, ...)
```

## Arguments

- x:

  An object containing a predictor definition

- ...:

  Arguments passed on to submethods.

## Value

`TRUE` if the expression is believed to be rowwise, `FALSE` otherwise.

## Examples

``` r
bru_is_rowwise(new_bru_pred_expr(~ x + y))
#> [1] FALSE
```
