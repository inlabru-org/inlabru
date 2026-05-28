# Output size of a mapping

Implementations must return an integer denoting the mapper output
length. The default implementation returns `NROW(input)`. Mappers such
as `bm_multi` and `bm_collect`, that can accept
[`list()`](https://rdrr.io/r/base/list.html) inputs require their own
method implementations.

## Usage

``` r
ibm_n_output(mapper, input, state = NULL, inla_f = FALSE, ...)

# Default S3 method
ibm_n_output(mapper, input, state = NULL, inla_f = FALSE, ...)

# S3 method for class 'bm_taylor'
ibm_n_output(mapper, input, ...)

# S3 method for class 'bm_shift'
ibm_n_output(mapper, input, state = NULL, ..., n_state = NULL)

# S3 method for class 'bm_scale'
ibm_n_output(mapper, input, state = NULL, ..., n_state = NULL)

# S3 method for class 'bm_aggregate'
ibm_n_output(mapper, input = NULL, ...)

# S3 method for class 'bm_marginal'
ibm_n_output(mapper, input, state = NULL, ..., n_state = NULL)

# S3 method for class 'bm_pipe'
ibm_n_output(mapper, input, state = NULL, ..., n_state = NULL)

# S3 method for class 'bm_multi'
ibm_n_output(mapper, input, ...)

# S3 method for class 'bm_reparam'
ibm_n_output(mapper, ...)

# S3 method for class 'bm_collect'
ibm_n_output(mapper, input, state = NULL, inla_f = FALSE, multi = FALSE, ...)

# S3 method for class 'bm_repeat'
ibm_n_output(mapper, ...)

# S3 method for class 'bm_sum'
ibm_n_output(mapper, input, state = NULL, ...)
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

- n_state:

  integer giving the length of the state vector for mappers that have
  state dependent output size.

- multi:

  logical; If `TRUE` (or positive), recurse one level into sub-mappers

## Methods (by class)

- `ibm_n_output(default)`: Returns `NROW(input)`

## See also

Other mapper methods:
[`bru_mapper_generics`](https://inlabru-org.github.io/inlabru/reference/bru_mapper_generics.md),
[`ibm_eval()`](https://inlabru-org.github.io/inlabru/reference/ibm_eval.md),
[`ibm_eval2()`](https://inlabru-org.github.io/inlabru/reference/ibm_eval2.md),
[`ibm_inla_subset()`](https://inlabru-org.github.io/inlabru/reference/ibm_inla_subset.md),
[`ibm_invalid_output()`](https://inlabru-org.github.io/inlabru/reference/ibm_invalid_output.md),
[`ibm_is_linear()`](https://inlabru-org.github.io/inlabru/reference/ibm_is_linear.md),
[`ibm_is_rowwise()`](https://inlabru-org.github.io/inlabru/reference/ibm_is_rowwise.md),
[`ibm_jacobian()`](https://inlabru-org.github.io/inlabru/reference/ibm_jacobian.md),
[`ibm_linear()`](https://inlabru-org.github.io/inlabru/reference/ibm_linear.md),
[`ibm_n()`](https://inlabru-org.github.io/inlabru/reference/ibm_n.md),
[`ibm_names()`](https://inlabru-org.github.io/inlabru/reference/ibm_names.md),
[`ibm_simplify()`](https://inlabru-org.github.io/inlabru/reference/ibm_simplify.md),
[`ibm_values()`](https://inlabru-org.github.io/inlabru/reference/ibm_values.md)
