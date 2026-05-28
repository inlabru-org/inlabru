# Find index subset of INLA visible states

Implementations must return a logical vector of `TRUE/FALSE` for the
subset such that, given the full A matrix and values output,
`A[, subset, drop = FALSE]` and `values[subset]` (or
`values[subset, , drop = FALSE]` for data.frame values) are equal to the
`inla_f = TRUE` version of A and values. The default method uses the
`ibm_values` output to construct the subset indexing.

## Usage

``` r
ibm_inla_subset(mapper, ...)

# Default S3 method
ibm_inla_subset(mapper, ...)
```

## Arguments

- mapper:

  A mapper S3 object, inheriting from `bru_mapper`.

- ...:

  Arguments passed on to other methods

## Methods (by class)

- `ibm_inla_subset(default)`: Uses the
  \[ibm_values`()] output to construct the inla subset indexing as the difference between `inla_f=FALSE`and`inla_f=TRUE`. Extra arguments such as `multi`are passed on to [ibm_values()]. This means it supports both regular vector values and`multi=1\`
  data.frame values.

## See also

Other mapper methods:
[`bru_mapper_generics`](https://inlabru-org.github.io/inlabru/reference/bru_mapper_generics.md),
[`ibm_eval()`](https://inlabru-org.github.io/inlabru/reference/ibm_eval.md),
[`ibm_eval2()`](https://inlabru-org.github.io/inlabru/reference/ibm_eval2.md),
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
