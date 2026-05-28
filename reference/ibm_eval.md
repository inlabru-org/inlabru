# Evaluate a mapping

Implementations must return a vector of length
[`ibm_n_output()`](https://inlabru-org.github.io/inlabru/reference/ibm_n_output.md).
The `input` contents must be in a format accepted by
[`ibm_jacobian()`](https://inlabru-org.github.io/inlabru/reference/ibm_jacobian.md)
for the mapper.

## Usage

``` r
ibm_eval(mapper, input, state = NULL, ...)

# Default S3 method
ibm_eval(mapper, input, state = NULL, ..., jacobian = NULL)

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

- jacobian:

  For `ibm_eval()` methods, an optional pre-computed Jacobian, typically
  supplied by internal methods that already have the Jacobian.

- log:

  logical; control `log` output. Default `TRUE`, see the `ibm_eval()`
  details for `logsumexp` mappers.

- logit:

  logical; control `logit` output. Default `TRUE`, see the `ibm_eval()`
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

- pre_A:

  **\[deprecated\]** in favour of `jacobian`.

- multi:

  logical; If `TRUE` (or positive), recurse one level into sub-mappers

- sub_lin:

  Internal, optional pre-computed sub-mapper information

## Methods (by class)

- `ibm_eval(default)`: Verifies that the mapper is linear with
  [`ibm_is_linear()`](https://inlabru-org.github.io/inlabru/reference/ibm_is_linear.md),
  and then computes a linear mapping as `ibm_jacobian(...) %*% state`.
  When `state` is `NULL`, a zero vector of length
  [`ibm_n_output()`](https://inlabru-org.github.io/inlabru/reference/ibm_n_output.md)
  is returned.

- `ibm_eval(bm_taylor)`: Evaluates linearised mapper information at the
  given `state`. The `input` argument is ignored, so that the usual
  argument order `ibm_eval(mapper, input, state)` syntax can be used,
  but also `ibm_eval(mapper, state = state)`. For a mapper with a named
  jacobian list, the `state` argument must also be a named list. If
  `state` is `NULL`, all-zero is assumed.

- `ibm_eval(bm_const)`: Returns the input values, with NA replaced by 0.

- `ibm_eval(bm_logsumexp)`: When `log` is `TRUE` (default), `ibm_eval()`
  for `logsumexp` returns the log-sum-weight-exp value. If `FALSE`, the
  `sum-weight-exp` value is returned.

- `ibm_eval(bm_logitaverage)`: When `logit` is `TRUE` (default),
  `ibm_eval()` for `logitaverage` returns the
  logit-sum-weight-inverse-logit value. If `FALSE`, the
  `sum-weights=invere-logit` value is returned.

- `ibm_eval(bm_marginal)`: When `xor(mapper[["inverse"]], reverse)` is
  `FALSE`, `ibm_eval()` for `marginal` returns `qfun(pnorm(x), param)`,
  evaluated in a numerically stable way. Otherwise, evaluates the
  inverse `qnorm(pfun(x, param))` instead.

## See also

Other mapper methods:
[`bru_mapper_generics`](https://inlabru-org.github.io/inlabru/reference/bru_mapper_generics.md),
[`ibm_eval2()`](https://inlabru-org.github.io/inlabru/reference/ibm_eval2.md),
[`ibm_inla_subset()`](https://inlabru-org.github.io/inlabru/reference/ibm_inla_subset.md),
[`ibm_invalid_output()`](https://inlabru-org.github.io/inlabru/reference/ibm_invalid_output.md),
[`ibm_is_linear()`](https://inlabru-org.github.io/inlabru/reference/ibm_is_linear.md),
[`ibm_is_rowwise()`](https://inlabru-org.github.io/inlabru/reference/ibm_is_rowwise.md),
[`ibm_jacobian()`](https://inlabru-org.github.io/inlabru/reference/ibm_jacobian.md),
[`ibm_linear()`](https://inlabru-org.github.io/inlabru/reference/ibm_linear.md),
[`ibm_n()`](https://inlabru-org.github.io/inlabru/reference/ibm_n.md),
[`ibm_n_output()`](https://inlabru-org.github.io/inlabru/reference/ibm_n_output.md),
[`ibm_names()`](https://inlabru-org.github.io/inlabru/reference/ibm_names.md),
[`ibm_simplify()`](https://inlabru-org.github.io/inlabru/reference/ibm_simplify.md),
[`ibm_values()`](https://inlabru-org.github.io/inlabru/reference/ibm_values.md)
