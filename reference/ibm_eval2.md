# Evaluate a mapper and its Jacobian

Implementations must return a list with elements `offset` and
`jacobian`. The `input` contents must be in a format accepted by
[`ibm_jacobian()`](https://inlabru-org.github.io/inlabru/reference/ibm_jacobian.md)
for the mapper.

## Usage

``` r
ibm_eval2(mapper, input, state = NULL, ...)

# Default S3 method
ibm_eval2(mapper, input, state = NULL, ...)

# S3 method for class 'bm_pipe'
ibm_eval2(mapper, input, state = NULL, ...)

# S3 method for class 'bm_expr'
ibm_eval2(mapper, input, state = NULL, ..., data = NULL)

# S3 method for class 'bru_obs'
ibm_eval2(
  mapper,
  input,
  state,
  ...,
  multi = FALSE,
  comp_mappers,
  eval_fun = NULL
)

# S3 method for class 'bru_obs_list'
ibm_eval2(
  mapper,
  input,
  state,
  ...,
  multi = FALSE,
  comp_mappers,
  eval_fun = NULL
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

A list with elements `offset` and `jacobian`, where `offset` is a vector
of length `ibm_n_output(mapper, input, state, ...)`, and `jacobian` is a
matrix of size `ibm_n_output(mapper, input, state, ...)` by
`ibm_n(mapper, inla_f = FALSE)`.

## Methods (by class)

- `ibm_eval2(default)`: Calls `jacobian <- ibm_jacobian(...)` and
  `offset <- ibm_eval(..., jacobian = jacobian)` and returns a list with
  elements `offset` and `jacobian`, as needed by
  [`ibm_as_taylor.default()`](https://inlabru-org.github.io/inlabru/reference/ibm_as_taylor.md)
  and similar methods. Mapper classes can implement their own
  `ibm_eval2` method if joint construction of evaluation and Jacobian is
  more efficient than separate or sequential construction.

## See also

Other mapper methods:
[`bru_mapper_generics`](https://inlabru-org.github.io/inlabru/reference/bru_mapper_generics.md),
[`ibm_as_taylor()`](https://inlabru-org.github.io/inlabru/reference/ibm_as_taylor.md),
[`ibm_eval()`](https://inlabru-org.github.io/inlabru/reference/ibm_eval.md),
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
ibm_eval2(m, input = c(1, 3, 4, 5, 2), state = 2)
#> $offset
#> [1]  2  6  8 10  4
#> 
#> $jacobian
#> 5 x 1 sparse Matrix of class "dgCMatrix"
#>       
#> [1,] 1
#> [2,] 3
#> [3,] 4
#> [4,] 5
#> [5,] 2
#> 
```
