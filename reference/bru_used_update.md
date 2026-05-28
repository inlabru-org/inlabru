# Update used_component information objects

Merge available component labels information with used components
information.

## Usage

``` r
bru_used_update(x, labels, ...)

# S3 method for class 'bru_obs_list'
bru_used_update(x, labels, ...)

# S3 method for class 'bru_obs'
bru_used_update(x, labels, ...)

# S3 method for class 'bru_pred_expr'
bru_used_update(x, labels, ...)

# S3 method for class 'bru_used'
bru_used_update(x, labels, ...)
```

## Arguments

- x:

  Object to be updated

- labels:

  character vector of component labels

- ...:

  Unused

## Value

An updated version of `x`

## See also

Other bru_used:
[`bru_used()`](https://inlabru-org.github.io/inlabru/reference/bru_used.md),
[`bru_used_vars()`](https://inlabru-org.github.io/inlabru/reference/bru_used_vars.md),
[`new_bru_used()`](https://inlabru-org.github.io/inlabru/reference/new_bru_used.md)
