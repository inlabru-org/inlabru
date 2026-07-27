# Iterated INLA

This is an internal wrapper for iterated runs of
[`INLA::inla`](https://rdrr.io/pkg/INLA/man/inla.html). For nonlinear
models, a linearisation is done with `bru_compute_linearisation`, with a
line search method between each iteration. The
[`INLA::inla.stack`](https://rdrr.io/pkg/INLA/man/inla.stack.html)
information is setup by
[`bru_make_stack()`](https://inlabru-org.github.io/inlabru/reference/bru_make_stack.md).

## Usage

``` r
iinla(
  model,
  initial = NULL,
  options,
  lhoods = deprecated(),
  inputs = deprecated()
)
```

## Arguments

- model:

  A
  [bru_model](https://inlabru-org.github.io/inlabru/reference/bru_model.md)
  object

- initial:

  A previous `bru` result or a list of named latent variable initial
  states (missing elements are set to zero), to be used as starting
  point, or `NULL`. If non-null, overrides `options$bru_initial`

- options:

  A `bru_options` object.

- lhoods:

  **\[deprecated\]** Deprecated from version `2.14.1.9011`, since the
  [bru_obs_list](https://inlabru-org.github.io/inlabru/reference/bru_obs.md)
  information is now part of the
  [bru_model](https://inlabru-org.github.io/inlabru/reference/bru_model.md)
  object.

- inputs:

  **\[deprecated\]** Deprecated from version `2.14.1.9011`, since the
  inputs information is now part of the
  [bru_model](https://inlabru-org.github.io/inlabru/reference/bru_model.md)
  object. Optional pre-computed list of per-likelihood component
  evaluations, from
  [`bru_input.bru_obs_list()`](https://inlabru-org.github.io/inlabru/reference/bru_input.md).

## Value

An `iinla` object that inherits from
[`INLA::inla`](https://rdrr.io/pkg/INLA/man/inla.html), with an added
field `bru_iinla` with elements

- log:

  The diagnostic log messages produced by the run

- states:

  The list of linearisation points, one for each inla run

- inla_stack:

  The `inla.stack` object from the final inla run

- track:

  A list of convergence tracking vectors

If an inla run is aborted by an error, the returned object also contains
an element `error` with the error object.
