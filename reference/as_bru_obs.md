# Conversion methods for `bru_obs` and `bru_obs_list` objects

Methods for converting to `bru_obs` and `bru_obs_list` objects.

## Usage

``` r
as_bru_obs(x, ...)

as_bru_obs_list(x, .tag = NULL)

# S3 method for class 'bru_obs'
as_bru_obs(x, ...)

# S3 method for class 'bru_obs'
as_bru_obs_list(x, .tag = NULL)

# S3 method for class 'list'
as_bru_obs_list(x, .tag = NULL)

# S3 method for class 'bru_obs_list'
as_bru_obs_list(x, .tag = NULL)

# S3 method for class 'bru'
as_bru_obs_list(x, .tag = NULL)

# S3 method for class 'bru_info'
as_bru_obs_list(x, .tag = NULL)

# S3 method for class 'bru_model'
as_bru_obs_list(x, .tag = NULL)
```

## Arguments

- x:

  An object to convert to
  [bru_obs](https://inlabru-org.github.io/inlabru/reference/bru_obs.md)
  or
  [bru_obs_list](https://inlabru-org.github.io/inlabru/reference/bru_obs.md)

- ...:

  Additional arguments passed to sub-methods.

- .tag:

  character; optional name for the single observation model in the
  returned
  [bru_obs_list](https://inlabru-org.github.io/inlabru/reference/bru_obs.md)
  object. Default is `NULL`, which results in automatic naming based on
  the `tag` attribute of `x`, if present.

## Value

An object of class
[bru_obs](https://inlabru-org.github.io/inlabru/reference/bru_obs.md) or
[bru_obs_list](https://inlabru-org.github.io/inlabru/reference/bru_obs.md).

## See also

[`as_bru_comp_list()`](https://inlabru-org.github.io/inlabru/reference/as_bru_comp.md)
