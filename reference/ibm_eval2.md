# Evaluate a mapper and its Jacobian

Implementations must return a list with elements `offset` and
`jacobian`. The `input` contents must be in a format accepted by
[`ibm_jacobian()`](https://inlabru-org.github.io/inlabru/reference/ibm_jacobian.md)
for the mapper.

## Usage

``` r
ibm_eval2(mapper, input, state = NULL, ...)

# Default S3 method
ibm_eval2(mapper, input, state, ...)

# S3 method for class 'bm_pipe'
ibm_eval2(mapper, input, state = NULL, ...)
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

## Methods (by class)

- `ibm_eval2(default)`: Calls `jacobian <- ibm_jacobian(...)` and
  `offset <- ibm_eval(..., jacobian = jacobian)` and returns a list with
  elements `offset` and `jacobian`, as needed by
  [`ibm_linear.default()`](https://inlabru-org.github.io/inlabru/reference/ibm_linear.md)
  and similar methods. Mapper classes can implement their own
  `ibm_eval2` method if joint construction of evaluation and Jacobian is
  more efficient than separate or sequential construction.

## See also

Other mapper methods:
[`bru_mapper_generics`](https://inlabru-org.github.io/inlabru/reference/bru_mapper_generics.md),
[`ibm_eval()`](https://inlabru-org.github.io/inlabru/reference/ibm_eval.md),
[`ibm_inla_subset()`](https://inlabru-org.github.io/inlabru/reference/ibm_inla_subset.md),
[`ibm_invalid_output()`](https://inlabru-org.github.io/inlabru/reference/ibm_invalid_output.md),
[`ibm_is_linear()`](https://inlabru-org.github.io/inlabru/reference/ibm_is_linear.md),
[`ibm_is_rowwise()`](https://inlabru-org.github.io/inlabru/reference/ibm_is_rowwise.md),
[`ibm_jacobian()`](https://inlabru-org.github.io/inlabru/reference/ibm_jacobian.md),
[`ibm_linear()`](https://inlabru-org.github.io/inlabru/reference/ibm_linear.md),
[`ibm_n()`](https://inlabru-org.github.io/inlabru/reference/ibm_n.md),
[`ibm_n_output()`](https://inlabru-org.github.io/inlabru/reference/ibm_n_output.md),
[`ibm_names()`](https://inlabru-org.github.io/inlabru/reference/ibm_names.md),
[`ibm_simplify()`](https://inlabru-org.github.io/inlabru/reference/ibm_simplify.md),
[`ibm_values()`](https://inlabru-org.github.io/inlabru/reference/ibm_values.md)
