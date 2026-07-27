# Mapper for tensor product domains

Constructs a row-wise Kronecker product mapping of linear/affine
mappers. Any offset in sub-mappers is added into a combined offset. Only
linear/affine sub-mappers are allowed.

## Usage

``` r
bm_multi(mappers, simplify = FALSE)

bru_mapper_multi(...)

# S3 method for class 'bm_multi'
x[i, drop = TRUE]

# S3 method for class 'bru_mapper_multi'
x[i, drop = TRUE]
```

## Arguments

- mappers:

  A list of `bru_mapper` objects

- simplify:

  logical; If `TRUE`, removes trivial submappers. Currently only
  sub-mappers of class
  [`bm_index()`](https://inlabru-org.github.io/inlabru/reference/bm_index.md)
  with `ibm_n() == 1L` are removed, and only if the mappers are named
  (to avoid ordering mismatches). Default: `FALSE`

- ...:

  Arguments passed on to `bm_multi()`

- x:

  object from which to extract element(s)

- i:

  indices specifying element(s) to extract

- drop:

  logical; For `[.bm_multi`, whether to extract an individual mapper
  when `i` identifies a single element. If `FALSE`, a list of
  sub-mappers is returned (suitable e.g. for creating a new `bm_multi`
  object). Default: `TRUE`

## Value

A `bm_multi` mapper object.

- `[`-indexing a `bm_multi` extracts a subset `bm_multi` object (for
  drop `FALSE`) or an individual sub-mapper (for drop `TRUE`, and `i`
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
[`bm_pipe()`](https://inlabru-org.github.io/inlabru/reference/bm_pipe.md),
[`bm_reparam()`](https://inlabru-org.github.io/inlabru/reference/bm_reparam.md),
[`bm_repeat()`](https://inlabru-org.github.io/inlabru/reference/bm_repeat.md),
[`bm_scale()`](https://inlabru-org.github.io/inlabru/reference/bm_scale.md),
[`bm_shift()`](https://inlabru-org.github.io/inlabru/reference/bm_shift.md),
[`bm_sum()`](https://inlabru-org.github.io/inlabru/reference/bm_sum.md),
[`bm_taylor()`](https://inlabru-org.github.io/inlabru/reference/bm_taylor.md),
[`bru_get_mapper()`](https://inlabru-org.github.io/inlabru/reference/bru_get_mapper.md),
[`bru_mapper()`](https://inlabru-org.github.io/inlabru/reference/bru_mapper.md)

## Examples

``` r
(m <- bm_multi(list(a = bm_index(2), b = bm_index(3))))
#> multi(a = {index}, b = {index})
ibm_eval2(m, list(a = c(1, 2, 1), b = c(1, 3, 2)), 1:6)
#> $offset
#> [1] 1 6 3
#> 
#> $jacobian
#> 3 x 6 sparse Matrix of class "dgCMatrix"
#>                 
#> [1,] 1 . . . . .
#> [2,] . . . . . 1
#> [3,] . . 1 . . .
#> 
```
