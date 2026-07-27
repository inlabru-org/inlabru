# Names of submapper

Implementations must return a character vector of sub-mapper names, or
`NULL`. Intended for providing information about multi-mappers and
mapper collections.

## Usage

``` r
ibm_names(mapper)

ibm_names(mapper) <- value

# Default S3 method
ibm_names(mapper, ...)

# S3 method for class 'bm_multi'
ibm_names(mapper)

# S3 method for class 'bm_multi'
ibm_names(mapper) <- value

# S3 method for class 'bru_mapper_multi'
ibm_names(mapper) <- value

# S3 method for class 'bm_collect'
ibm_names(mapper)

# S3 method for class 'bm_collect'
ibm_names(mapper) <- value

# S3 method for class 'bru_mapper_collect'
ibm_names(mapper) <- value

# S3 method for class 'bm_sum'
ibm_names(mapper)

# S3 method for class 'bm_sum'
ibm_names(mapper) <- value

# S3 method for class 'bru_mapper_sum'
ibm_names(mapper) <- value
```

## Arguments

- mapper:

  A mapper S3 object, inheriting from `bru_mapper`.

- value:

  a character vector of up to the same length as the number of mappers
  in the multi-mapper x

- ...:

  Arguments passed on to other methods

## Value

A character vector or `NULL`

## Methods (by class)

- `ibm_names(default)`: Returns `NULL`

- `ibm_names(bm_multi)`: Returns the names from the sub-mappers list

- `ibm_names(bm_collect)`: Returns the names from the sub-mappers list

- `ibm_names(bm_sum)`: Returns the names from the sub-mappers list

## Functions

- `ibm_names(mapper) <- value`: Set mapper names.

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
[`ibm_simplify()`](https://inlabru-org.github.io/inlabru/reference/ibm_simplify.md),
[`ibm_values()`](https://inlabru-org.github.io/inlabru/reference/ibm_values.md)

## Examples

``` r
# ibm_names
mapper <- bm_multi(list(A = bm_index(2), B = bm_index(2)))
ibm_names(mapper)
#> [1] "A" "B"
ibm_names(mapper) <- c("new", "names")
ibm_names(mapper)
#> [1] "new"   "names"
m <- bm_multi(list(A = bm_linear(), B = bm_linear()))
ibm_names(m)
#> [1] "A" "B"
```
