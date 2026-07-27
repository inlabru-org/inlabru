# Simplify a mapper

Implementations must return a
[bru_mapper](https://inlabru-org.github.io/inlabru/reference/bru_mapper.md)
object. The default method returns the result of
[`ibm_as_taylor()`](https://inlabru-org.github.io/inlabru/reference/ibm_as_taylor.md)
for linear/affine mappers, and the original `mapper` for non-linear
mappers.

## Usage

``` r
ibm_simplify(mapper, input = NULL, state = NULL, ...)

# Default S3 method
ibm_simplify(mapper, input = NULL, state = NULL, ...)

# S3 method for class 'bm_pipe'
ibm_simplify(
  mapper,
  input = NULL,
  state = NULL,
  inla_f = FALSE,
  ...,
  n_state = NULL
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

- n_state:

  integer giving the length of the state vector for mappers that have
  state dependent output size.

## Value

The original mapper is returned for non-linear mappers, and the output
of
[`ibm_as_taylor()`](https://inlabru-org.github.io/inlabru/reference/ibm_as_taylor.md)
is returned for linear mappers.

## Methods (by class)

- `ibm_simplify(default)`: Calls
  [`ibm_as_taylor()`](https://inlabru-org.github.io/inlabru/reference/ibm_as_taylor.md)
  for linear mappers, and returns the original mapper for non-linear
  mappers.

- `ibm_simplify(bm_pipe)`: Constructs a simplified `pipe` mapper. For
  fully linear pipes, calls
  [`ibm_as_taylor()`](https://inlabru-org.github.io/inlabru/reference/ibm_as_taylor.md).
  For partially non-linear pipes, replaces each sequence of linear
  mappers with a single
  [`bm_taylor()`](https://inlabru-org.github.io/inlabru/reference/bm_taylor.md)
  mapper, while keeping the full list of original mapper names, allowing
  the original `input` structure to be used also with the simplified
  mappers, since the `taylor` mappers are not dependent on inputs.

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
[`ibm_values()`](https://inlabru-org.github.io/inlabru/reference/ibm_values.md)

## Examples

``` r
m <- bm_linear()
ibm_simplify(m, input = c(1, 3, 4, 5, 2), state = 2)
#> taylor
```
