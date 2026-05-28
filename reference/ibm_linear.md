# Compute a mapper linearisation

Implementations must return a
[bm_taylor](https://inlabru-org.github.io/inlabru/reference/bm_taylor.md)
object The linearisation information includes `offset`, `jacobian`, and
`state0`. The state information indicates for which state the `offset`
was evaluated, with `NULL` meaning all-zero. The linearised mapper
output is defined as

    effect(input, state) =
      offset(input, state0) + jacobian(input, state0) %*% (state - state0)

The default method calls
[`ibm_eval()`](https://inlabru-org.github.io/inlabru/reference/ibm_eval.md)
and
[`ibm_jacobian()`](https://inlabru-org.github.io/inlabru/reference/ibm_jacobian.md)
to generate the needed information.

## Usage

``` r
ibm_linear(mapper, input, state = NULL, ...)

# Default S3 method
ibm_linear(mapper, input, state, ...)

# S3 method for class 'bm_multi'
ibm_linear(mapper, input, state, inla_f = FALSE, ...)

# S3 method for class 'bm_collect'
ibm_linear(mapper, input, state, inla_f = FALSE, ...)

# S3 method for class 'bm_repeat'
ibm_linear(mapper, input, state, ...)

# S3 method for class 'bm_sum'
ibm_linear(mapper, input, state, ...)

# S3 method for class 'bru_comp'
ibm_linear(mapper, input, state = NULL, ...)
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

## Value

A
[bm_taylor](https://inlabru-org.github.io/inlabru/reference/bm_taylor.md)
object. The `state0` information in the affine mapper indicates for
which state the `offset` was evaluated; The affine mapper output is
defined as

    effect(input, state) =
      offset(input, state0) + jacobian(input, state0) %*% (state - state0)

## Methods (by class)

- `ibm_linear(default)`: Calls
  [`ibm_eval2()`](https://inlabru-org.github.io/inlabru/reference/ibm_eval2.md)
  and returns a
  [bm_taylor](https://inlabru-org.github.io/inlabru/reference/bm_taylor.md)
  object.

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
[`ibm_n()`](https://inlabru-org.github.io/inlabru/reference/ibm_n.md),
[`ibm_n_output()`](https://inlabru-org.github.io/inlabru/reference/ibm_n_output.md),
[`ibm_names()`](https://inlabru-org.github.io/inlabru/reference/ibm_names.md),
[`ibm_simplify()`](https://inlabru-org.github.io/inlabru/reference/ibm_simplify.md),
[`ibm_values()`](https://inlabru-org.github.io/inlabru/reference/ibm_values.md)
