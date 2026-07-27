# Evaluate a mapping

Implementations must return a vector of length
[`ibm_n_output()`](https://inlabru-org.github.io/inlabru/reference/ibm_n_output.md).
The `input` contents must be in a format accepted by
[`ibm_jacobian()`](https://inlabru-org.github.io/inlabru/reference/ibm_jacobian.md)
for the mapper.

Specific implementations for
[`bm_aggregate`](https://inlabru-org.github.io/inlabru/reference/ibm_eval_methods.md),
[`bm_collect`](https://inlabru-org.github.io/inlabru/reference/ibm_eval_methods.md),
[`bm_const`](https://inlabru-org.github.io/inlabru/reference/ibm_eval_methods.md),
[`bm_expr`](https://inlabru-org.github.io/inlabru/reference/ibm_eval_methods.md),
[`bm_list`](https://inlabru-org.github.io/inlabru/reference/bru_model_mapper_methods.md),
[`bm_logitaverage`](https://inlabru-org.github.io/inlabru/reference/ibm_eval_methods.md),
[`bm_logsumexp`](https://inlabru-org.github.io/inlabru/reference/ibm_eval_methods.md),
[`bm_marginal`](https://inlabru-org.github.io/inlabru/reference/ibm_eval_methods.md),
[`bm_multi`](https://inlabru-org.github.io/inlabru/reference/ibm_eval_methods.md),
[`bm_pipe`](https://inlabru-org.github.io/inlabru/reference/ibm_eval_methods.md),
[`bm_reparam`](https://inlabru-org.github.io/inlabru/reference/ibm_eval_methods.md),
[`bm_repeat`](https://inlabru-org.github.io/inlabru/reference/ibm_eval_methods.md),
[`bm_scale`](https://inlabru-org.github.io/inlabru/reference/ibm_eval_methods.md),
[`bm_shift`](https://inlabru-org.github.io/inlabru/reference/ibm_eval_methods.md),
[`bm_sum`](https://inlabru-org.github.io/inlabru/reference/ibm_eval_methods.md),
[`bm_taylor`](https://inlabru-org.github.io/inlabru/reference/ibm_eval_methods.md),
[`bru_comp`](https://inlabru-org.github.io/inlabru/reference/bru_model_mapper_methods.md),
[`bru_comp_list`](https://inlabru-org.github.io/inlabru/reference/bru_model_mapper_methods.md),
[`bru_model`](https://inlabru-org.github.io/inlabru/reference/bru_model_mapper_methods.md),
`bru_obs`, `bru_obs_list`, `default`.

## Usage

``` r
ibm_eval(mapper, input, state = NULL, ...)

# Default S3 method
ibm_eval(mapper, input, state = NULL, ..., jacobian = NULL)

# S3 method for class 'bru_obs'
ibm_eval(
  mapper,
  input,
  state,
  ...,
  multi = FALSE,
  comp_mappers,
  eval_fun = NULL
)

# S3 method for class 'bru_obs_list'
ibm_eval(
  mapper,
  input,
  state,
  ...,
  multi = FALSE,
  comp_mappers,
  eval_fun = NULL
)
```

## Arguments

- mapper:

  A mapper S3 object, inheriting from `bru_mapper`.

- input:

  Data input for the mapper.

- state:

  A vector of latent state values for the mapping, of length
  `ibm_n(mapper, inla_f = FALSE)`

- ...:

  Arguments passed on to
  [`ibm_eval_methods`](https://inlabru-org.github.io/inlabru/reference/ibm_eval_methods.md)

  `log`

  :   logical; control `log` output. Default `TRUE`, see the
      `ibm_eval()` details for `logsumexp` mappers.

  `logit`

  :   logical; control `logit` output. Default `TRUE`, see the
      `ibm_eval()` details for `logitaverage` mappers.

  `reverse`

  :   logical; control `bm_marginal` evaluation. Default `FALSE`. When
      `TRUE`, reverses the direction of the mapping, see details for
      `marginal` mappers.

  `pre_A`

  :   **\[deprecated\]** in favour of `jacobian`.

  `sub_lin`

  :   Internal, optional pre-computed sub-mapper information

  `data_mask`

  :   A data mask object to use for evaluating the expression. If `NULL`
      or missing, a data mask will be constructed with
      [`bru_data_mask()`](https://inlabru-org.github.io/inlabru/reference/bru_data_mask.md)
      from the `data`, `state`, `derived`, and `.envir` arguments. This
      can be used to avoid redundant construction of the data mask when
      evaluating different expressions multiple times with the same
      input and state.

  `.envir`

  :   The environment in which to evaluate the expression. By default,
      this is set to the caller environment.

  `derived`

  :   The state vectors of variables derived from the root variables. If
      `NULL` or missing, defaults to an empty list.

  `inla_f`

  :   logical; when `TRUE` for
      [`ibm_n()`](https://inlabru-org.github.io/inlabru/reference/ibm_n.md)
      and
      [`ibm_values()`](https://inlabru-org.github.io/inlabru/reference/ibm_values.md),
      the result must be compatible with the `INLA::f(...)` and
      corresponding `INLA::inla.stack(...)` constructions. For
      `ibm_{eval,jacobian,linear}`, the `input` interpretation may be
      different. Implementations do not normally need to do anything
      different, except for mappers of the type needed for hidden
      multicomponent models such as "bym2", which can be handled by
      `bm_collect`.

  `data`

  :   should be a list with data objects, with the main object called
      `data`; see
      [`bm_expr()`](https://inlabru-org.github.io/inlabru/reference/bm_expr.md)
      for details.

- jacobian:

  For `ibm_eval()` methods, an optional pre-computed Jacobian, typically
  supplied by internal methods that already have the Jacobian.

- multi:

  logical; If `TRUE` (or positive), recurse one level into sub-mappers

- comp_mappers:

  A list of mappers, typically from `as_bm_list<bru_comp_list>`.

- eval_fun:

  A list of functions, typically from
  [`bru_eval_fun()`](https://inlabru-org.github.io/inlabru/reference/bru_eval_fun.md).

## Value

A vector of length `ibm_n_output(mapper, input, state, ...)`.

## Methods (by class)

- `ibm_eval(default)`: Verifies that the mapper is linear with
  [`ibm_is_linear()`](https://inlabru-org.github.io/inlabru/reference/ibm_is_linear.md),
  and then computes a linear mapping as `ibm_jacobian(...) %*% state`.
  When `state` is `NULL`, a zero vector of length
  [`ibm_n_output()`](https://inlabru-org.github.io/inlabru/reference/ibm_n_output.md)
  is returned.

## See also

Other mapper methods:
[`bru_mapper_generics`](https://inlabru-org.github.io/inlabru/reference/bru_mapper_generics.md),
[`ibm_as_taylor()`](https://inlabru-org.github.io/inlabru/reference/ibm_as_taylor.md),
[`ibm_eval2()`](https://inlabru-org.github.io/inlabru/reference/ibm_eval2.md),
[`ibm_inla_subset()`](https://inlabru-org.github.io/inlabru/reference/ibm_inla_subset.md),
[`ibm_invalid_output()`](https://inlabru-org.github.io/inlabru/reference/ibm_invalid_output.md),
[`ibm_is_linear()`](https://inlabru-org.github.io/inlabru/reference/ibm_is_linear.md),
[`ibm_is_rowwise()`](https://inlabru-org.github.io/inlabru/reference/ibm_is_rowwise.md),
[`ibm_jacobian()`](https://inlabru-org.github.io/inlabru/reference/ibm_jacobian.md),
[`ibm_n()`](https://inlabru-org.github.io/inlabru/reference/ibm_n.md),
[`ibm_n_output()`](https://inlabru-org.github.io/inlabru/reference/ibm_n_output.md),
[`ibm_names()`](https://inlabru-org.github.io/inlabru/reference/ibm_names.md),
[`ibm_simplify()`](https://inlabru-org.github.io/inlabru/reference/ibm_simplify.md),
[`ibm_values()`](https://inlabru-org.github.io/inlabru/reference/ibm_values.md)

## Examples

``` r
m <- bm_linear()
ibm_eval(m, input = c(1, 3, 4, 5, 2), state = 2)
#> [1]  2  6  8 10  4
```
