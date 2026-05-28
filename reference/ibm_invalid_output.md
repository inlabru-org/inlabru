# Detect invalid input to a mapper

Implementations should return a logical vector of length
`ibm_n_output(mapper, input, state, ...)` indicating which, if any,
output elements of `ibm_eval(mapper, input, state, ...)` are known to be
invalid. For for multi/collect mappers, a list, when given a
`multi=TRUE` argument.

## Usage

``` r
ibm_invalid_output(mapper, input, state, ...)

# Default S3 method
ibm_invalid_output(mapper, input, state, ...)

# S3 method for class 'bm_index'
ibm_invalid_output(mapper, input, state, ...)

# S3 method for class 'bm_multi'
ibm_invalid_output(mapper, input, state, inla_f = FALSE, multi = FALSE, ...)

# S3 method for class 'bm_collect'
ibm_invalid_output(mapper, input, state, inla_f = FALSE, multi = FALSE, ...)

# S3 method for class 'bm_repeat'
ibm_invalid_output(mapper, input, state, ...)

# S3 method for class 'bm_sum'
ibm_invalid_output(mapper, input, state, multi = FALSE, ...)
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

- multi:

  logical; If `TRUE` (or positive), recurse one level into sub-mappers

## Methods (by class)

- `ibm_invalid_output(default)`: Returns an all-`FALSE` logical vector.

- `ibm_invalid_output(bm_index)`: Returns `TRUE` out-of-range input
  values. Returns `FALSE` for in-range and `NA` input values (to support
  NA giving zero effect).

- `ibm_invalid_output(bm_multi)`: Accepts a list with named entries, or
  a list with unnamed but ordered elements. The names must match the
  sub-mappers, see
  [`ibm_names.bm_multi()`](https://inlabru-org.github.io/inlabru/reference/ibm_names.md).
  Each list element should take a format accepted by the corresponding
  sub-mapper. In case each element is a vector, the input can be given
  as a data.frame with named columns, a matrix with named columns, or a
  matrix with unnamed but ordered columns.

- `ibm_invalid_output(bm_collect)`: Accepts a list with named entries,
  or a list with unnamed but ordered elements. The names must match the
  sub-mappers, see
  [`ibm_names.bm_collect()`](https://inlabru-org.github.io/inlabru/reference/ibm_names.md).
  Each list element should take a format accepted by the corresponding
  sub-mapper. In case each element is a vector, the input can be given
  as a data.frame with named columns, a matrix with named columns, or a
  matrix with unnamed but ordered columns.

- `ibm_invalid_output(bm_repeat)`: Passes on the input to the
  corresponding method.

- `ibm_invalid_output(bm_sum)`: Accepts a list with named entries, or a
  list with unnamed but ordered elements. The names must match the
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
[`ibm_is_linear()`](https://inlabru-org.github.io/inlabru/reference/ibm_is_linear.md),
[`ibm_is_rowwise()`](https://inlabru-org.github.io/inlabru/reference/ibm_is_rowwise.md),
[`ibm_jacobian()`](https://inlabru-org.github.io/inlabru/reference/ibm_jacobian.md),
[`ibm_linear()`](https://inlabru-org.github.io/inlabru/reference/ibm_linear.md),
[`ibm_n()`](https://inlabru-org.github.io/inlabru/reference/ibm_n.md),
[`ibm_n_output()`](https://inlabru-org.github.io/inlabru/reference/ibm_n_output.md),
[`ibm_names()`](https://inlabru-org.github.io/inlabru/reference/ibm_names.md),
[`ibm_simplify()`](https://inlabru-org.github.io/inlabru/reference/ibm_simplify.md),
[`ibm_values()`](https://inlabru-org.github.io/inlabru/reference/ibm_values.md)
