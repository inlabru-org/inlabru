# Value vector for a mapping

When `inla_f=TRUE`, implementations must return a vector that would be
interpretable by an `INLA::f(..., values = ...)` specification. The
exception is the method for `bm_multi`, that returns a multi-column data
frame if `multi=TRUE`.

## Usage

``` r
ibm_values(mapper, inla_f = FALSE, ...)

# Default S3 method
ibm_values(mapper, inla_f = FALSE, ...)

# S3 method for class 'bm_fmesher'
ibm_values(mapper, ...)

# S3 method for class 'bm_fm_mesh_1d'
ibm_values(mapper, ...)

# S3 method for class 'bm_taylor'
ibm_values(mapper, inla_f = FALSE, multi = FALSE, ...)

# S3 method for class 'bm_linear'
ibm_values(mapper, ...)

# S3 method for class 'bm_matrix'
ibm_values(mapper, ...)

# S3 method for class 'bm_factor'
ibm_values(mapper, ...)

# S3 method for class 'bm_const'
ibm_values(mapper, ...)

# S3 method for class 'bm_shift'
ibm_values(mapper, ..., state = NULL, n_state = NULL)

# S3 method for class 'bm_scale'
ibm_values(mapper, ..., state = NULL, n_state = NULL)

# S3 method for class 'bm_aggregate'
ibm_values(mapper, ..., state = NULL, n_state = NULL)

# S3 method for class 'bm_marginal'
ibm_values(mapper, ..., state = NULL, n_state = NULL)

# S3 method for class 'bm_pipe'
ibm_values(mapper, ...)

# S3 method for class 'bm_multi'
ibm_values(mapper, inla_f = FALSE, multi = FALSE, ...)

# S3 method for class 'bm_reparam'
ibm_values(mapper, ...)

# S3 method for class 'bm_collect'
ibm_values(mapper, inla_f = FALSE, multi = FALSE, ...)

# S3 method for class 'bm_expr'
ibm_values(mapper, inla_f = FALSE, ...)

# S3 method for class 'bm_repeat'
ibm_values(mapper, ...)

# S3 method for class 'bm_sum'
ibm_values(mapper, inla_f = FALSE, multi = FALSE, ...)
```

## Arguments

- mapper:

  A mapper S3 object, inheriting from `bru_mapper`.

- inla_f:

  logical; when `TRUE` for
  [`ibm_n()`](https://inlabru-org.github.io/inlabru/reference/ibm_n.md)
  and `ibm_values()`, the result must be compatible with the
  `INLA::f(...)` and corresponding `INLA::inla.stack(...)`
  constructions. For `ibm_{eval,jacobian,linear}`, the `input`
  interpretation may be different. Implementations do not normally need
  to do anything different, except for mappers of the type needed for
  hidden multicomponent models such as "bym2", which can be handled by
  `bm_collect`.

- ...:

  Arguments passed on to other methods

- multi:

  logical; If `TRUE` (or positive), recurse one level into sub-mappers

- state:

  A vector of latent state values for the mapping, of length
  `ibm_n(mapper, inla_f = FALSE)`

- n_state:

  integer giving the length of the state vector for mappers that have
  state dependent output size.

## Value

A vector of length `ibm_n(mapper, inla_f = FALSE)`

## Methods (by class)

- `ibm_values(default)`: Returns a non-null element 'values' from the
  mapper object, and `seq_len(ibm_n(mapper))` if it doesn't exist.

- `ibm_values(bm_fmesher)`: Returns an index vector for the mesh basis
  functions.

- `ibm_values(bm_fm_mesh_1d)`: Returns an index vector into the basis
  functions for an indexed mapper. Otherwise, the `mid` values if
  present in the mesh being mapped, and otherwise returns the `loc`
  values of the mesh.

- `ibm_values(bm_linear)`: Returns `1.0`

- `ibm_values(bm_matrix)`: For integer labels, the vector of labels. For
  character labels, the labels as a factor variable.

- `ibm_values(bm_factor)`: Returns the factor levels (minus the first
  level for `factor_mapping` "contrast"), or an integer vector (if
  `indexed = TRUE` in
  [`bm_factor()`](https://inlabru-org.github.io/inlabru/reference/bm_factor.md)).

## See also

Other mapper methods:
[`bru_mapper_generics`](https://inlabru-org.github.io/inlabru/reference/bru_mapper_generics.md),
[`ibm_as_taylor()`](https://inlabru-org.github.io/inlabru/reference/ibm_as_taylor.md),
[`ibm_eval()`](https://inlabru-org.github.io/inlabru/reference/ibm_eval.md),
[`ibm_eval2()`](https://inlabru-org.github.io/inlabru/reference/ibm_eval2.md),
[`ibm_inla_subset()`](https://inlabru-org.github.io/inlabru/reference/ibm_inla_subset.md),
[`ibm_invalid_output()`](https://inlabru-org.github.io/inlabru/reference/ibm_invalid_output.md),
[`ibm_is_linear()`](https://inlabru-org.github.io/inlabru/reference/ibm_is_linear.md),
[`ibm_is_rowwise()`](https://inlabru-org.github.io/inlabru/reference/ibm_is_rowwise.md),
[`ibm_jacobian()`](https://inlabru-org.github.io/inlabru/reference/ibm_jacobian.md),
[`ibm_n()`](https://inlabru-org.github.io/inlabru/reference/ibm_n.md),
[`ibm_n_output()`](https://inlabru-org.github.io/inlabru/reference/ibm_n_output.md),
[`ibm_names()`](https://inlabru-org.github.io/inlabru/reference/ibm_names.md),
[`ibm_simplify()`](https://inlabru-org.github.io/inlabru/reference/ibm_simplify.md)

## Examples

``` r
m <- bm_index(4)
ibm_values(m)
#> [1] 1 2 3 4
```
