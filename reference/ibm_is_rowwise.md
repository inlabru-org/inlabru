# Check if a mapper is rowwise

Implementations must return `TRUE` or `FALSE`. If `TRUE` (returned by
the default method unless the mapper contains an `is_rowwise` variable),
users of the mapper may assume the mapper uses its inputs in "rowwise"
manner, so that blockwise evaluation is always possible.

## Usage

``` r
ibm_is_rowwise(mapper, ...)

# Default S3 method
ibm_is_rowwise(mapper, ...)
```

## Arguments

- mapper:

  A mapper S3 object, inheriting from `bru_mapper`.

- ...:

  Arguments passed on to other methods

## Methods (by class)

- `ibm_is_rowwise(default)`: Returns logical `is_rowwise` from the
  mapper object if it exists, and otherwise `TRUE`.

## See also

Other mapper methods:
[`bru_mapper_generics`](https://inlabru-org.github.io/inlabru/reference/bru_mapper_generics.md),
[`ibm_eval()`](https://inlabru-org.github.io/inlabru/reference/ibm_eval.md),
[`ibm_eval2()`](https://inlabru-org.github.io/inlabru/reference/ibm_eval2.md),
[`ibm_inla_subset()`](https://inlabru-org.github.io/inlabru/reference/ibm_inla_subset.md),
[`ibm_invalid_output()`](https://inlabru-org.github.io/inlabru/reference/ibm_invalid_output.md),
[`ibm_is_linear()`](https://inlabru-org.github.io/inlabru/reference/ibm_is_linear.md),
[`ibm_jacobian()`](https://inlabru-org.github.io/inlabru/reference/ibm_jacobian.md),
[`ibm_linear()`](https://inlabru-org.github.io/inlabru/reference/ibm_linear.md),
[`ibm_n()`](https://inlabru-org.github.io/inlabru/reference/ibm_n.md),
[`ibm_n_output()`](https://inlabru-org.github.io/inlabru/reference/ibm_n_output.md),
[`ibm_names()`](https://inlabru-org.github.io/inlabru/reference/ibm_names.md),
[`ibm_simplify()`](https://inlabru-org.github.io/inlabru/reference/ibm_simplify.md),
[`ibm_values()`](https://inlabru-org.github.io/inlabru/reference/ibm_values.md)
