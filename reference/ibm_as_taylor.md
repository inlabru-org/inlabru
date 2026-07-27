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
[`ibm_eval2()`](https://inlabru-org.github.io/inlabru/reference/ibm_eval2.md)
to generate the needed information.

## Usage

``` r
ibm_as_taylor(mapper, input, state = NULL, ...)

# Default S3 method
ibm_as_taylor(mapper, input, state, ...)

# S3 method for class 'bm_multi'
ibm_as_taylor(mapper, input, state, inla_f = FALSE, ...)

# S3 method for class 'bm_collect'
ibm_as_taylor(mapper, input, state, inla_f = FALSE, ...)

# S3 method for class 'bm_expr'
ibm_as_taylor(
  mapper,
  input,
  state,
  ...,
  derived = NULL,
  jacobians = NULL,
  data = data,
  inla_f = FALSE
)

# S3 method for class 'bru_obs'
ibm_as_taylor(mapper, input, state, ..., multi = FALSE, comp_mappers)

# S3 method for class 'bru_obs_list'
ibm_as_taylor(
  mapper,
  input,
  state,
  ...,
  multi = FALSE,
  comp_mappers,
  eval_fun = NULL
)

# S3 method for class 'bm_repeat'
ibm_as_taylor(mapper, input, state, ...)

# S3 method for class 'bm_sum'
ibm_as_taylor(mapper, input, state, ...)

# S3 method for class 'bru_comp'
ibm_as_taylor(mapper, input, state = NULL, ...)
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

- derived:

  The state vectors of variables derived from the root variables. If
  `NULL` or missing, defaults to an empty list.

- jacobians:

  The state vectors of variables derived from the root If `jacobians` is
  `NULL` or missing, defaults to an empty list. If not `NULL`, should be
  a list with named entries, one for variable derived from the root
  variables. Each list element should be a named list of Jacobian
  matrices, with names matching the root variables. Missing entries are
  treated as all-zero matrices.

- data:

  should be a list with data objects, with the main object called
  `data`; see
  [`bm_expr()`](https://inlabru-org.github.io/inlabru/reference/bm_expr.md)
  for details.

- multi:

  logical; If `TRUE` (or positive), recurse one level into sub-mappers

- comp_mappers:

  A list of mappers, typically from `as_bm_list<bru_comp_list>`.

- eval_fun:

  A list of functions, typically from
  [`bru_eval_fun()`](https://inlabru-org.github.io/inlabru/reference/bru_eval_fun.md).

## Value

A
[bm_taylor](https://inlabru-org.github.io/inlabru/reference/bm_taylor.md)
object. The `state0` information in the affine mapper indicates for
which state the `offset` was evaluated; The affine mapper output is
defined as

    effect(input, state) =
      offset(input, state0) + jacobian(input, state0) %*% (state - state0)

## Methods (by class)

- `ibm_as_taylor(default)`: Calls
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

## Examples

``` r
m <- bm_linear()
ibm_as_taylor(m, input = c(1, 3, 4, 5, 2), state = 2)
#> taylor
```
