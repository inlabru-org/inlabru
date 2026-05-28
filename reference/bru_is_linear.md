# Check for predictor linearity

Checks if a predictor expression is linear (or affine)

## Usage

``` r
bru_is_linear(x, ...)

# S3 method for class 'bru'
bru_is_linear(x, ...)

# S3 method for class 'bru_info'
bru_is_linear(x, ...)

# S3 method for class 'bru_model'
bru_is_linear(x, ...)

# S3 method for class 'bru_obs'
bru_is_linear(x, ...)

# S3 method for class 'bru_obs_list'
bru_is_linear(x, ...)

# S3 method for class 'bru_pred_expr'
bru_is_linear(x, ...)

# S3 method for class 'bru_comp_list'
bru_is_linear(x, ...)

# S3 method for class 'bru_comp'
bru_is_linear(x, ...)

# S3 method for class 'bru_mapper'
bru_is_linear(x, ...)
```

## Arguments

- x:

  A object containing a predictor definition

- ...:

  Arguments passed on recursively.

## Value

`TRUE` if the expression is detected to be linear, `FALSE` otherwise.

## Examples

``` r
bru_is_linear(new_bru_pred_expr(~ x + y))
#> [1] TRUE
bru_is_linear(new_bru_pred_expr(~ x * y))
#> [1] FALSE
bru_is_linear(bm_scale())
#> [1] TRUE
bru_is_linear(bm_logsumexp())
#> [1] FALSE
```
