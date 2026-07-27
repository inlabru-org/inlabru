# Mapper for adding multiple mappers

Defines a mapper that adds the effects of each submapper. The
[`ibm_n()`](https://inlabru-org.github.io/inlabru/reference/ibm_n.md)
method returns the sum of `ibm_n(mappers[[k]])`, and
[`ibm_values()`](https://inlabru-org.github.io/inlabru/reference/ibm_values.md)
returns `seq_len(ibm_n(mapper))`.

## Usage

``` r
bm_sum(mappers, single_input = FALSE)

bru_mapper_sum(...)

# S3 method for class 'bm_sum'
x[i, drop = TRUE]

# S3 method for class 'bru_mapper_sum'
x[i, drop = TRUE]
```

## Arguments

- mappers:

  A list of `bru_mapper` objects.

- single_input:

  logical. If `TRUE`, the input is passed to all sub-mappers. Otherwise,
  the input should be a list, data.frame, or matrix. If the `mappers`
  list has named entries, the `input` can reference their corresponding
  sub-mapper using its name.

- ...:

  Arguments passed on to `bm_sum()`

- x:

  object from which to extract element(s)

- i:

  indices specifying element(s) to extract

- drop:

  logical; For `[.bm_sum`, whether to extract an individual mapper when
  `i` identifies a single element. If `FALSE`, a list of sub-mappers is
  returned (suitable e.g. for creating a new `bm_sum` object). Default:
  `TRUE`

## Value

A `bm_sum` object.

- `[`-indexing a `bm_sum` extracts a subset `bm_sum` object (for drop
  `FALSE`) or an individual sub-mapper (for drop `TRUE`, and `i`
  identifies a single element)

## See also

[bru_mapper](https://inlabru-org.github.io/inlabru/reference/bru_mapper.md),
[bru_mapper_generics](https://inlabru-org.github.io/inlabru/reference/bru_mapper_generics.md)

Other mappers:
[`bm_aggregate()`](https://inlabru-org.github.io/inlabru/reference/bm_aggregate.md),
[`bm_collect()`](https://inlabru-org.github.io/inlabru/reference/bm_collect.md),
[`bm_const()`](https://inlabru-org.github.io/inlabru/reference/bm_const.md),
[`bm_expr()`](https://inlabru-org.github.io/inlabru/reference/bm_expr.md),
[`bm_factor()`](https://inlabru-org.github.io/inlabru/reference/bm_factor.md),
[`bm_fm_mesh_1d`](https://inlabru-org.github.io/inlabru/reference/bm_fm_mesh_1d.md),
[`bm_fmesher()`](https://inlabru-org.github.io/inlabru/reference/bm_fmesher.md),
[`bm_harmonics()`](https://inlabru-org.github.io/inlabru/reference/bm_harmonics.md),
[`bm_index()`](https://inlabru-org.github.io/inlabru/reference/bm_index.md),
[`bm_linear()`](https://inlabru-org.github.io/inlabru/reference/bm_linear.md),
[`bm_logitaverage()`](https://inlabru-org.github.io/inlabru/reference/bm_logitaverage.md),
[`bm_logsumexp()`](https://inlabru-org.github.io/inlabru/reference/bm_logsumexp.md),
[`bm_marginal()`](https://inlabru-org.github.io/inlabru/reference/bm_marginal.md),
[`bm_matrix()`](https://inlabru-org.github.io/inlabru/reference/bm_matrix.md),
[`bm_multi()`](https://inlabru-org.github.io/inlabru/reference/bm_multi.md),
[`bm_pipe()`](https://inlabru-org.github.io/inlabru/reference/bm_pipe.md),
[`bm_reparam()`](https://inlabru-org.github.io/inlabru/reference/bm_reparam.md),
[`bm_repeat()`](https://inlabru-org.github.io/inlabru/reference/bm_repeat.md),
[`bm_scale()`](https://inlabru-org.github.io/inlabru/reference/bm_scale.md),
[`bm_shift()`](https://inlabru-org.github.io/inlabru/reference/bm_shift.md),
[`bm_taylor()`](https://inlabru-org.github.io/inlabru/reference/bm_taylor.md),
[`bru_get_mapper()`](https://inlabru-org.github.io/inlabru/reference/bru_get_mapper.md),
[`bru_mapper()`](https://inlabru-org.github.io/inlabru/reference/bru_mapper.md)

## Examples

``` r
(m <- bm_sum(list(a = bm_index(3), b = bm_index(2))))
#> sum(a = {index}, b = {index})
ibm_n(m)
#> [1] 5
ibm_values(m)
#> [1] 1 2 3 4 5
ibm_jacobian(m, list(a = 1:3, b = c(1, 1, 2)))
#> 3 x 5 sparse Matrix of class "dgCMatrix"
#>               
#> [1,] 1 . . 1 .
#> [2,] . 1 . 1 .
#> [3,] . . 1 . 1
ibm_eval(
  m,
  list(a = 1:3, b = c(1, 1, 2)),
  seq_len(ibm_n(m))
)
#> [1] 5 6 8
```
