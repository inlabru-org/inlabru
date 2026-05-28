# Utility functions for bru observation model objects

Utility functions for bru observation model objects

## Usage

``` r
bru_obs_inla_family(x, ...)

# S3 method for class 'bru_obs'
bru_obs_inla_family(x, ...)

# S3 method for class 'bru_obs_list'
bru_obs_inla_family(x, ...)

bru_obs_control_family(x, control.family = NULL, ...)

# S3 method for class 'bru_obs'
bru_obs_control_family(x, control.family = NULL, ...)

# S3 method for class 'bru_obs_list'
bru_obs_control_family(x, control.family = NULL, ...)

bru_obs_control_gcpo(x, ...)

# S3 method for class 'bru_obs'
bru_obs_control_gcpo(
  x,
  index_offset,
  index_length = bru_response_size(x),
  force_weights,
  ...
)

# S3 method for class 'bru_obs_list'
bru_obs_control_gcpo(x, control.gcpo = NULL, ...)
```

## Arguments

- x:

  Object of `bru_obs` or `bru_obs_list` type

- ...:

  Further arguments passed on to the submethods

- control.family:

  list of INLA `control.family` options to override

- index_offset:

  integer; offset to add to indices in `control.gcpo`

- index_length:

  integer; length of the response vector for the observation model.
  Defaults to `bru_response_size(x)`.

- force_weights:

  logical; if `TRUE`, ensure that `control.gcpo$weights` is populated.
  This is needed if any of the observation models in a `bru_obs_list`
  has non-null `control.gcpo$weights`.

- control.gcpo:

  list of INLA `control.gcpo` default options

## Value

- `bru_obs_inla_family()` returns a string or vector of strings

&nbsp;

- `bru_obs_control_family()` returns a list with
  [`INLA::control.family`](https://rdrr.io/pkg/INLA/man/control.family.html)
  options, or a list of such lists, with one element per observation
  model

&nbsp;

- `bru_obs_control_gcpo()` returns a list with
  [`INLA::control.gcpo`](https://rdrr.io/pkg/INLA/man/control.gcpo.html)
  options, with predictor/response variable indices unified for
  multi-observation models. If a `bru_obs` model has a NULL
  `control.gcpo` argument, an empty list is returned for that model.

## See also

[`summary.bru_obs()`](https://inlabru-org.github.io/inlabru/reference/bru_obs_print.md)
