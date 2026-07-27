# Size of the latent vector of a mapping

Implementations must return the size of the latent vector being mapped.

## Usage

``` r
ibm_n(mapper, inla_f = FALSE, ...)

# Default S3 method
ibm_n(mapper, inla_f = FALSE, ...)

# S3 method for class 'bm_fmesher'
ibm_n(mapper, ...)

# S3 method for class 'bm_fm_mesh_1d'
ibm_n(mapper, ...)

# S3 method for class 'bm_taylor'
ibm_n(mapper, inla_f = FALSE, multi = FALSE, ...)

# S3 method for class 'bm_linear'
ibm_n(mapper, ...)

# S3 method for class 'bm_matrix'
ibm_n(mapper, ...)

# S3 method for class 'bm_factor'
ibm_n(mapper, ...)

# S3 method for class 'bm_const'
ibm_n(mapper, ...)

# S3 method for class 'bm_shift'
ibm_n(mapper, ..., state = NULL, n_state = NULL)

# S3 method for class 'bm_scale'
ibm_n(mapper, ..., state = NULL, n_state = NULL)

# S3 method for class 'bm_aggregate'
ibm_n(mapper, ..., input = NULL, state = NULL, n_state = NULL)

# S3 method for class 'bm_marginal'
ibm_n(mapper, ..., state = NULL, n_state = NULL)

# S3 method for class 'bm_pipe'
ibm_n(mapper, ..., input = NULL, state = NULL)

# S3 method for class 'bm_multi'
ibm_n(mapper, inla_f = FALSE, multi = FALSE, ...)

# S3 method for class 'bm_harmonics'
ibm_n(mapper, inla_f = FALSE, ...)

# S3 method for class 'bm_reparam'
ibm_n(mapper, ...)

# S3 method for class 'bm_collect'
ibm_n(mapper, inla_f = FALSE, multi = FALSE, ...)

# S3 method for class 'bm_expr'
ibm_n(mapper, ..., input = NULL, state = NULL, multi = FALSE, data = NULL)

# S3 method for class 'bm_repeat'
ibm_n(mapper, ...)

# S3 method for class 'bm_sum'
ibm_n(mapper, inla_f = FALSE, multi = FALSE, ...)
```

## Arguments

- mapper:

  A mapper S3 object, inheriting from `bru_mapper`.

- inla_f:

  logical; when `TRUE` for `ibm_n()` and
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

- state:

  A vector of latent state values for the mapping, of length
  `ibm_n(mapper, inla_f = FALSE)`

- n_state:

  integer giving the length of the state vector for mappers that have
  state dependent output size.

- input:

  Data input for the mapper.

- data:

  should be a list with data objects, with the main object called
  `data`; see
  [`bm_expr()`](https://inlabru-org.github.io/inlabru/reference/bm_expr.md)
  for details.

## Value

An integer denoting the size of the latent vector being mapped.

## Methods (by class)

- `ibm_n(default)`: Returns a non-null element 'n' from the mapper
  object, and gives an error if it doesn't exist. If `inla_f=TRUE`,
  first checks for a 'n_inla' element.

- `ibm_n(bm_fmesher)`: Returns the
  [`fmesher::fm_dof()`](https://inlabru-org.github.io/fmesher/reference/fm_dof.html)
  value of the mesh being mapped.

- `ibm_n(bm_fm_mesh_1d)`: Returns the
  [`fmesher::fm_dof()`](https://inlabru-org.github.io/fmesher/reference/fm_dof.html)
  value of the mesh being mapped.

- `ibm_n(bm_linear)`: Returns `1L`

- `ibm_n(bm_matrix)`: Returns the number of columns in the matrix
  mapper.

- `ibm_n(bm_factor)`: Returns the number of levels when `factor_mapping`
  is "full", and the number of levels minus one if `factor_mapping` is
  "contrast".

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
[`ibm_n_output()`](https://inlabru-org.github.io/inlabru/reference/ibm_n_output.md),
[`ibm_names()`](https://inlabru-org.github.io/inlabru/reference/ibm_names.md),
[`ibm_simplify()`](https://inlabru-org.github.io/inlabru/reference/ibm_simplify.md),
[`ibm_values()`](https://inlabru-org.github.io/inlabru/reference/ibm_values.md)
