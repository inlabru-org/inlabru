# Utility functions for bru observation model GCPO

Extract and combine
[`INLA::control.gcpo`](https://rdrr.io/pkg/INLA/man/control.gcpo.html)
options from `bru_obs` and `bru_obs_list` objects, adjusting observation
indices for multi-likelihood models where the response vector is the
concatenation of all likelihoods.

## Usage

``` r
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

  A `bru_obs` or `bru_obs_list` object.

- ...:

  Further arguments passed to submethods.

- index_offset:

  integer; offset to add to indices in `control.gcpo`, equal to the
  total number of response rows from all preceding likelihoods.

- index_length:

  integer; number of response rows for this likelihood. Defaults to
  `bru_response_size(x)`.

- force_weights:

  logical; if `TRUE`, ensure that `control.gcpo$weights` is populated.
  Required when any likelihood in a `bru_obs_list` has non-NULL
  `control.gcpo$weights`.

- control.gcpo:

  list of
  [`INLA::control.gcpo`](https://rdrr.io/pkg/INLA/man/control.gcpo.html)
  default options to merge with per-likelihood settings.

## Value

- `bru_obs_control_gcpo()` returns a list with
  [`INLA::control.gcpo`](https://rdrr.io/pkg/INLA/man/control.gcpo.html)
  options, with predictor/response variable indices unified for
  multi-observation models. If a `bru_obs` model has a
  `control.gcpo = NULL` argument, an empty list is returned for that
  model.

## See also

[`bru_obs()`](https://inlabru-org.github.io/inlabru/reference/bru_obs.md),
[`bru_block_gcpo()`](https://inlabru-org.github.io/inlabru/reference/bru_block_gcpo.md)
