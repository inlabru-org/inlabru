# Compute inlabru model linearisation information

Compute inlabru model linearisation information

## Usage

``` r
bru_compute_linearisation(...)

# S3 method for class 'bru_comp'
bru_compute_linearisation(
  cmp,
  model,
  lhood_expr,
  data,
  data_extra,
  input,
  state,
  comp_simple,
  effects,
  pred0,
  used,
  is_rowwise,
  eps,
  n_pred = NULL,
  ...
)

# S3 method for class 'bru_obs'
bru_compute_linearisation(
  lhood,
  model,
  data,
  input,
  state,
  comp_simple,
  eps,
  ...
)

# S3 method for class 'bru_obs_list'
bru_compute_linearisation(
  lhoods,
  model,
  input,
  state,
  comp_simple,
  eps = 1e-05,
  ...
)

# S3 method for class 'bru_model'
bru_compute_linearisation(model, lhoods, input, state, comp_simple, ...)
```

## Arguments

- ...:

  Parameters passed on to other methods

- cmp:

  A
  [bru_comp](https://inlabru-org.github.io/inlabru/reference/bru_comp.md)
  object

- model:

  A `bru_model` object

- lhood_expr:

  A predictor expression

- data:

  Input data

- data_extra:

  Additional data for the predictor

- input:

  Precomputed component inputs from
  [`bru_input()`](https://inlabru-org.github.io/inlabru/reference/bru_input.md)

- state:

  The state information, as a list of named vectors

- comp_simple:

  Component evaluation information

  - For `bru_comp`: A
    [bm_taylor](https://inlabru-org.github.io/inlabru/reference/bm_taylor.md)
    object

  - For `bru_obs`: A
    [bm_list](https://inlabru-org.github.io/inlabru/reference/bm_list.md)
    object for the components in the likelihood

  - For `bru_obs_list`: A list of
    [bm_list](https://inlabru-org.github.io/inlabru/reference/bm_list.md)
    objects

- effects:

  - For `bru_comp`: Precomputed effect list for all components involved
    in the likelihood expression

- pred0:

  Precomputed predictor for the given state

- used:

  A
  [`bru_used()`](https://inlabru-org.github.io/inlabru/reference/bru_used.md)
  object for the predictor expression

- is_rowwise:

  logical; If `FALSE`, the predictor expression may involve several rows
  of the input data to influence the same row.

- eps:

  The finite difference step size

- n_pred:

  The length of the predictor expression. If not `NULL`, scalar
  predictor evaluations are expanded to vectors of length `n_pred`.

- lhood:

  A `bru_obs` object

- lhoods:

  A `bru_obs_list` object
