# Conversion methods for `bru_comp` and `bru_comp_list` objects

Methods for converting to `bru_comp` and `bru_comp_list` objects.

## Usage

``` r
as_bru_comp(x, ...)

as_bru_comp_list(x, ...)

# S3 method for class 'bru_comp'
as_bru_comp(x, ...)

# S3 method for class 'bru_comp_list'
as_bru_comp_list(x, ...)

# S3 method for class 'bru_comp'
as_bru_comp_list(x, ...)

# S3 method for class 'bru'
as_bru_comp_list(x, ...)

# S3 method for class 'bru_info'
as_bru_comp_list(x, ...)

# S3 method for class 'bru_model'
as_bru_comp_list(x, ...)

# S3 method for class 'list'
as_bru_comp_list(x, ...)

# S3 method for class 'formula'
as_bru_comp_list(x, ...)
```

## Arguments

- x:

  An object to convert to
  [bru_comp](https://inlabru-org.github.io/inlabru/reference/bru_comp.md)
  or
  [bru_comp_list](https://inlabru-org.github.io/inlabru/reference/bru_comp_list.md)

- ...:

  Additional arguments passed on to
  [`bru_comp_list()`](https://inlabru-org.github.io/inlabru/reference/bru_comp_list.md).

## Value

An object of class
[bru_comp_list](https://inlabru-org.github.io/inlabru/reference/bru_comp_list.md).

## Functions

- `as_bru_comp_list(bru)`: Extract the component list from a
  [`bru()`](https://inlabru-org.github.io/inlabru/reference/bru.md)
  object.

- `as_bru_comp_list(bru_info)`: Extract the component list from a
  [`bru_info()`](https://inlabru-org.github.io/inlabru/reference/bru_info.md)
  object.

- `as_bru_comp_list(bru_model)`: Extract the component list from a
  [`bru_model()`](https://inlabru-org.github.io/inlabru/reference/bru_model.md)
  object.

## See also

[`as_bru_obs()`](https://inlabru-org.github.io/inlabru/reference/as_bru_obs.md),
[`as_bru_obs_list()`](https://inlabru-org.github.io/inlabru/reference/as_bru_obs.md)
