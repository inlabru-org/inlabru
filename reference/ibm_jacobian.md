# Jacobian of a mapper

Implementations must return a (sparse) matrix of size
`ibm_n_output(mapper, input, inla_f)` by
`ibm_n(mapper, inla_f = FALSE)`. The `inla_f=TRUE` argument should only
affect the allowed type of input format.

## Usage

``` r
ibm_jacobian(mapper, input, state = NULL, inla_f = FALSE, ...)

# Default S3 method
ibm_jacobian(mapper, input, state = NULL, ...)

# S3 method for class 'bm_fmesher'
ibm_jacobian(mapper, input, ...)

# S3 method for class 'bm_fm_mesh_1d'
ibm_jacobian(mapper, input, ...)

# S3 method for class 'bm_index'
ibm_jacobian(mapper, input, state, ...)

# S3 method for class 'bm_taylor'
ibm_jacobian(mapper, ..., multi = FALSE)

# S3 method for class 'bm_linear'
ibm_jacobian(mapper, input, ...)

# S3 method for class 'bm_matrix'
ibm_jacobian(mapper, input, state = NULL, inla_f = FALSE, ...)

# S3 method for class 'bm_factor'
ibm_jacobian(mapper, input, ...)

# S3 method for class 'bm_const'
ibm_jacobian(mapper, input, ...)

# S3 method for class 'bm_shift'
ibm_jacobian(mapper, input, state = NULL, ...)

# S3 method for class 'bm_scale'
ibm_jacobian(mapper, input, state = NULL, ...)

# S3 method for class 'bm_aggregate'
ibm_jacobian(mapper, input, state = NULL, ...)

# S3 method for class 'bm_logsumexp'
ibm_jacobian(mapper, input, state = NULL, ...)

# S3 method for class 'bm_logitaverage'
ibm_jacobian(mapper, input, state = NULL, ...)

# S3 method for class 'bm_marginal'
ibm_jacobian(mapper, input, state = NULL, ..., reverse = FALSE)

# S3 method for class 'bm_pipe'
ibm_jacobian(mapper, input, state = NULL, ...)

# S3 method for class 'bm_multi'
ibm_jacobian(
  mapper,
  input,
  state = NULL,
  inla_f = FALSE,
  multi = FALSE,
  ...,
  sub_A = NULL
)

# S3 method for class 'bm_harmonics'
ibm_jacobian(mapper, input, state = NULL, inla_f = FALSE, ...)

# S3 method for class 'bm_reparam'
ibm_jacobian(mapper, input, state = NULL, ...)

# S3 method for class 'bm_collect'
ibm_jacobian(
  mapper,
  input,
  state = NULL,
  inla_f = FALSE,
  multi = FALSE,
  ...,
  sub_lin = NULL
)

# S3 method for class 'bm_repeat'
ibm_jacobian(
  mapper,
  input,
  state = NULL,
  inla_f = FALSE,
  multi = FALSE,
  ...,
  sub_lin = NULL
)

# S3 method for class 'bm_sum'
ibm_jacobian(
  mapper,
  input,
  state = NULL,
  inla_f = FALSE,
  multi = FALSE,
  ...,
  sub_lin = NULL
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

- ...:

  Arguments passed on to other methods

- multi:

  logical; If `TRUE` (or positive), recurse one level into sub-mappers

- reverse:

  logical; control `bm_marginal` evaluation. Default `FALSE`. When
  `TRUE`, reverses the direction of the mapping, see details for
  `marginal` mappers.

- sub_A:

  Internal; precomputed Jacobian matrices.

- sub_lin:

  Internal, optional pre-computed sub-mapper information

## Methods (by class)

- `ibm_jacobian(default)`: Mapper classes must implement their own
  `ibm_jacobian` method.

- `ibm_jacobian(bm_fmesher)`: Returns the
  [`fmesher::fm_basis()`](https://inlabru-org.github.io/fmesher/reference/fm_basis.html)
  matrix of the mesh being mapped.

- `ibm_jacobian(bm_fm_mesh_1d)`: Returns the
  [`fmesher::fm_basis()`](https://inlabru-org.github.io/fmesher/reference/fm_basis.html)
  matrix of the mesh being mapped.

- `ibm_jacobian(bm_matrix)`: Accepts `input` as a matrix, `Matrix`,
  `Spatial`, or `sfc_POINT` object.

- `ibm_jacobian(bm_shift)`: `input` NULL values are interpreted as no
  shift.

- `ibm_jacobian(bm_scale)`: `input` NULL values are interpreted as no
  scaling.

- `ibm_jacobian(bm_aggregate)`: `input` should be a list with elements
  `block` and `weights`. `block` should be a vector of the same length
  as the `state`, or `NULL`, with `NULL` equivalent to all-1. If
  `weights` is `NULL`, it's interpreted as all-1.

- `ibm_jacobian(bm_logsumexp)`: `input` should be a list with elements
  `block` and `weights`. `block` should be a vector of the same length
  as the `state`, or `NULL`, with `NULL` equivalent to all-1. If
  `weights` is `NULL`, it's interpreted as all-1.

- `ibm_jacobian(bm_logitaverage)`: `input` should be a list with
  elements `block` and `weights`. `block` should be a vector of the same
  length as the `state`, or `NULL`, with `NULL` equivalent to all-1. If
  `weights` is `NULL`, it's interpreted as all-1.

- `ibm_jacobian(bm_marginal)`: Non-NULL `input` values are interpreted
  as a parameter list for `qfun`, overriding that of the mapper itself.

- `ibm_jacobian(bm_multi)`: Accepts a list with named entries, or a list
  with unnamed but ordered elements. The names must match the
  sub-mappers, see
  [`ibm_names.bm_multi()`](https://inlabru-org.github.io/inlabru/reference/ibm_names.md).
  Each list element should take a format accepted by the corresponding
  sub-mapper. In case each element is a vector, the input can be given
  as a data.frame with named columns, a matrix with named columns, or a
  matrix with unnamed but ordered columns.

- `ibm_jacobian(bm_collect)`: Accepts a list with named entries, or a
  list with unnamed but ordered elements. The names must match the
  sub-mappers, see
  [`ibm_names.bm_collect()`](https://inlabru-org.github.io/inlabru/reference/ibm_names.md).
  Each list element should take a format accepted by the corresponding
  sub-mapper. In case each element is a vector, the input can be given
  as a data.frame with named columns, a matrix with named columns, or a
  matrix with unnamed but ordered columns. When `inla_f=TRUE` and
  `hidden=TRUE` in the mapper definition, the input format should
  instead match that of the first, non-hidden, sub-mapper.

- `ibm_jacobian(bm_repeat)`: The input should take the format of the
  repeated submapper.

- `ibm_jacobian(bm_sum)`: Accepts a list with named entries, or a list
  with unnamed but ordered elements. The names must match the
  sub-mappers, see
  [`ibm_names.bm_sum()`](https://inlabru-org.github.io/inlabru/reference/ibm_names.md).
  Each list element should take a format accepted by the corresponding
  sub-mapper. In case each element is a vector, the input can be given
  as a data.frame with named columns, a matrix with named columns, or a
  matrix with unnamed but ordered columns.

## See also

Other mapper methods:
[`bru_mapper_generics`](https://inlabru-org.github.io/inlabru/reference/bru_mapper_generics.md),
[`ibm_eval()`](https://inlabru-org.github.io/inlabru/reference/ibm_eval.md),
[`ibm_eval2()`](https://inlabru-org.github.io/inlabru/reference/ibm_eval2.md),
[`ibm_inla_subset()`](https://inlabru-org.github.io/inlabru/reference/ibm_inla_subset.md),
[`ibm_invalid_output()`](https://inlabru-org.github.io/inlabru/reference/ibm_invalid_output.md),
[`ibm_is_linear()`](https://inlabru-org.github.io/inlabru/reference/ibm_is_linear.md),
[`ibm_is_rowwise()`](https://inlabru-org.github.io/inlabru/reference/ibm_is_rowwise.md),
[`ibm_linear()`](https://inlabru-org.github.io/inlabru/reference/ibm_linear.md),
[`ibm_n()`](https://inlabru-org.github.io/inlabru/reference/ibm_n.md),
[`ibm_n_output()`](https://inlabru-org.github.io/inlabru/reference/ibm_n_output.md),
[`ibm_names()`](https://inlabru-org.github.io/inlabru/reference/ibm_names.md),
[`ibm_simplify()`](https://inlabru-org.github.io/inlabru/reference/ibm_simplify.md),
[`ibm_values()`](https://inlabru-org.github.io/inlabru/reference/ibm_values.md)
