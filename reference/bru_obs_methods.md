# Utility functions for bru observation model objects

Utility functions for bru observation model objects

## Usage

``` r
bru_obs_inla_family(x, ...)

# S3 method for class 'bru_obs'
bru_obs_inla_family(x, ...)

# S3 method for class 'bru_obs_list'
bru_obs_inla_family(x, ...)

# S3 method for class 'bru'
bru_obs_inla_family(x, ...)

bru_obs_family(x, ...)

# S3 method for class 'bru_obs'
bru_obs_family(x, ...)

# S3 method for class 'bru_obs_list'
bru_obs_family(x, ...)

# S3 method for class 'bru'
bru_obs_family(x, ...)

bru_obs_control_family(x, control.family = NULL, ...)

# S3 method for class 'bru_obs'
bru_obs_control_family(x, control.family = NULL, ...)

# S3 method for class 'bru_obs_list'
bru_obs_control_family(x, control.family = NULL, ...)
```

## Arguments

- x:

  Object of `bru_obs` or `bru_obs_list` type

- ...:

  Further arguments passed on to the submethods

- control.family:

  list of INLA `control.family` options to override

## Value

- `bru_obs_inla_family()` returns a string or vector of strings of the
  `family` name(s) used in the
  [`INLA::inla()`](https://rdrr.io/pkg/INLA/man/inla.html) call for the
  observation model(s) in `x`.

&nbsp;

- `bru_obs_family()` returns a string or vector of strings of the
  `family` name(s) used to define each
  [`bru_obs()`](https://inlabru-org.github.io/inlabru/reference/bru_obs.md).
  This may be different the internal technical name(s) used in the
  [`INLA::inla()`](https://rdrr.io/pkg/INLA/man/inla.html) call.

&nbsp;

- `bru_obs_control_family()` returns a list with
  [`INLA::control.family`](https://rdrr.io/pkg/INLA/man/control.family.html)
  options, or a list of such lists, with one element per observation
  model

## See also

[`summary.bru_obs()`](https://inlabru-org.github.io/inlabru/reference/bru_obs_print.md)
