# Check if a mapper is linear/affine

Implementations must return `TRUE` or `FALSE`. If `TRUE` (returned by
the default method unless the mapper contains an `is_linear` variable),
users of the mapper may assume the mapper is linear/affine.

## Usage

``` r
ibm_is_linear(mapper, ...)

# Default S3 method
ibm_is_linear(mapper, ...)

# S3 method for class 'bm_multi'
ibm_is_linear(mapper, multi = FALSE, ...)

# S3 method for class 'bm_collect'
ibm_is_linear(mapper, inla_f = FALSE, multi = FALSE, ...)

# S3 method for class 'bm_expr'
ibm_is_linear(mapper, ...)

# S3 method for class 'bm_sum'
ibm_is_linear(mapper, multi = FALSE, ...)
```

## Arguments

- mapper:

  A mapper S3 object, inheriting from `bru_mapper`.

- ...:

  Arguments passed on to other methods

- multi:

  logical; If `TRUE` (or positive), recurse one level into sub-mappers

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

## Value

logical; `TRUE` if the mapper is linear/affine, `FALSE` otherwise.

## Methods (by class)

- `ibm_is_linear(default)`: Returns logical `is_linear` from the mapper
  object if it exists, and otherwise `TRUE`.

## See also

Other mapper methods:
[`bru_mapper_generics`](https://inlabru-org.github.io/inlabru/reference/bru_mapper_generics.md),
[`ibm_as_taylor()`](https://inlabru-org.github.io/inlabru/reference/ibm_as_taylor.md),
[`ibm_eval()`](https://inlabru-org.github.io/inlabru/reference/ibm_eval.md),
[`ibm_eval2()`](https://inlabru-org.github.io/inlabru/reference/ibm_eval2.md),
[`ibm_inla_subset()`](https://inlabru-org.github.io/inlabru/reference/ibm_inla_subset.md),
[`ibm_invalid_output()`](https://inlabru-org.github.io/inlabru/reference/ibm_invalid_output.md),
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
ibm_is_linear(m)
#> [1] TRUE
```
