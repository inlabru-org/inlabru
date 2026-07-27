# Specific `ibm_eval` method implementations

Specific
[`ibm_eval()`](https://inlabru-org.github.io/inlabru/reference/ibm_eval.md)
method implementations.

## Usage

``` r
# S3 method for class 'bm_taylor'
ibm_eval(mapper, input = NULL, state = NULL, ...)

# S3 method for class 'bm_const'
ibm_eval(mapper, input, state = NULL, ...)

# S3 method for class 'bm_shift'
ibm_eval(mapper, input, state = NULL, ...)

# S3 method for class 'bm_scale'
ibm_eval(mapper, input, state = NULL, ...)

# S3 method for class 'bm_aggregate'
ibm_eval(mapper, input, state = NULL, ...)

# S3 method for class 'bm_logsumexp'
ibm_eval(mapper, input, state = NULL, log = TRUE, ...)

# S3 method for class 'bm_logitaverage'
ibm_eval(mapper, input, state = NULL, logit = TRUE, ...)

# S3 method for class 'bm_marginal'
ibm_eval(mapper, input, state = NULL, ..., reverse = FALSE)

# S3 method for class 'bm_pipe'
ibm_eval(mapper, input, state = NULL, ...)

# S3 method for class 'bm_multi'
ibm_eval(
  mapper,
  input,
  state = NULL,
  inla_f = FALSE,
  ...,
  jacobian = NULL,
  pre_A = deprecated()
)

# S3 method for class 'bm_reparam'
ibm_eval(mapper, input, state = NULL, ..., jacobian = NULL)

# S3 method for class 'bm_collect'
ibm_eval(
  mapper,
  input,
  state,
  inla_f = FALSE,
  multi = FALSE,
  ...,
  sub_lin = NULL
)

# S3 method for class 'bm_expr'
ibm_eval(
  mapper,
  input,
  state = NULL,
  ...,
  derived = NULL,
  data = NULL,
  data_mask = NULL,
  .envir = rlang::caller_env()
)

# S3 method for class 'bm_repeat'
ibm_eval(mapper, input, state, multi = FALSE, ..., sub_lin = NULL)

# S3 method for class 'bm_sum'
ibm_eval(mapper, input, state, multi = FALSE, ..., sub_lin = NULL)
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

  Arguments passed on to other methods

- log:

  logical; control `log` output. Default `TRUE`, see the
  [`ibm_eval()`](https://inlabru-org.github.io/inlabru/reference/ibm_eval.md)
  details for `logsumexp` mappers.

- logit:

  logical; control `logit` output. Default `TRUE`, see the
  [`ibm_eval()`](https://inlabru-org.github.io/inlabru/reference/ibm_eval.md)
  details for `logitaverage` mappers.

- reverse:

  logical; control `bm_marginal` evaluation. Default `FALSE`. When
  `TRUE`, reverses the direction of the mapping, see details for
  `marginal` mappers.

- inla_f:

  logical; when `TRUE` for
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

- jacobian:

  For
  [`ibm_eval()`](https://inlabru-org.github.io/inlabru/reference/ibm_eval.md)
  methods, an optional pre-computed Jacobian, typically supplied by
  internal methods that already have the Jacobian.

- pre_A:

  **\[deprecated\]** in favour of `jacobian`.

- multi:

  logical; If `TRUE` (or positive), recurse one level into sub-mappers

- sub_lin:

  Internal, optional pre-computed sub-mapper information

- derived:

  The state vectors of variables derived from the root variables. If
  `NULL` or missing, defaults to an empty list.

- data:

  should be a list with data objects, with the main object called
  `data`; see
  [`bm_expr()`](https://inlabru-org.github.io/inlabru/reference/bm_expr.md)
  for details.

- data_mask:

  A data mask object to use for evaluating the expression. If `NULL` or
  missing, a data mask will be constructed with
  [`bru_data_mask()`](https://inlabru-org.github.io/inlabru/reference/bru_data_mask.md)
  from the `data`, `state`, `derived`, and `.envir` arguments. This can
  be used to avoid redundant construction of the data mask when
  evaluating different expressions multiple times with the same input
  and state.

- .envir:

  The environment in which to evaluate the expression. By default, this
  is set to the caller environment.

## Functions

- `ibm_eval(bm_taylor)`: Evaluates linearised mapper information at the
  given `state`. The `input` argument is ignored, so that the usual
  argument order `ibm_eval(mapper, input, state)` syntax can be used,
  but also `ibm_eval(mapper, state = state)`. For a mapper with a named
  jacobian list, the `state` argument must also be a named list. If
  `state` is `NULL`, all-zero is assumed.

- `ibm_eval(bm_const)`: Returns the input values, with NA replaced by 0.

- `ibm_eval(bm_logsumexp)`: When `log` is `TRUE` (default),
  [`ibm_eval()`](https://inlabru-org.github.io/inlabru/reference/ibm_eval.md)
  for `logsumexp` returns the log-sum-weight-exp value. If `FALSE`, the
  `sum-weight-exp` value is returned.

- `ibm_eval(bm_logitaverage)`: When `logit` is `TRUE` (default),
  [`ibm_eval()`](https://inlabru-org.github.io/inlabru/reference/ibm_eval.md)
  for `logitaverage` returns the logit-sum-weight-inverse-logit value.
  If `FALSE`, the `sum-weights=invere-logit` value is returned.

- `ibm_eval(bm_marginal)`: When `xor(mapper[["inverse"]], reverse)` is
  `FALSE`,
  [`ibm_eval()`](https://inlabru-org.github.io/inlabru/reference/ibm_eval.md)
  for `marginal` returns `qfun(pnorm(x), param)`, evaluated in a
  numerically stable way. Otherwise, evaluates the inverse
  `qnorm(pfun(x, param))` instead.

- `ibm_eval(bm_expr)`: Accepts a `state` list with named entries, one
  for each variable. The `input` and `data` formats should match the
  description given for
  [`bm_expr()`](https://inlabru-org.github.io/inlabru/reference/bm_expr.md).
