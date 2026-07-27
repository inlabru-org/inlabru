# Mapper for concatenated variables

Constructs a concatenated collection mapping

## Usage

``` r
bm_collect(mappers, hidden = FALSE)

bru_mapper_collect(...)

# S3 method for class 'bm_collect'
x[i, drop = TRUE]

# S3 method for class 'bru_mapper_collect'
x[i, drop = TRUE]
```

## Arguments

- mappers:

  A list of `bru_mapper` objects

- hidden:

  `logical`, set to `TRUE` to flag that the mapper is to be used as a
  first level input mapper for
  [`INLA::f()`](https://rdrr.io/pkg/INLA/man/f.html) in a model that
  requires making only the first mapper visible to
  [`INLA::f()`](https://rdrr.io/pkg/INLA/man/f.html) and
  [`INLA::inla.stack()`](https://rdrr.io/pkg/INLA/man/inla.stack.html),
  such as for "bym2" models, as activated by the `inla_f` argument to
  `ibm_n`, `ibm_values`, and `ibm_jacobian`. Set to `FALSE` to always
  access the full mapper, e.g. for `rgeneric` models

- ...:

  Arguments passed on to
  [`bm_scale()`](https://inlabru-org.github.io/inlabru/reference/bm_scale.md)

- x:

  object from which to extract element(s)

- i:

  indices specifying element(s) to extract

- drop:

  logical; For `[.bm_collect`, whether to extract an individual mapper
  when `i` identifies a single element. If `FALSE`, a list of
  sub-mappers is returned (suitable e.g. for creating a new `bm_collect`
  object). Default: `TRUE`

## Value

- `[`-indexing a `bm_collect` extracts a subset `bm_collect` object (for
  drop `FALSE`) or an individual sub-mapper (for drop `TRUE`, and `i`
  identifies a single element)

## See also

[bru_mapper](https://inlabru-org.github.io/inlabru/reference/bru_mapper.md),
[bru_mapper_generics](https://inlabru-org.github.io/inlabru/reference/bru_mapper_generics.md)

Other mappers:
[`bm_aggregate()`](https://inlabru-org.github.io/inlabru/reference/bm_aggregate.md),
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
[`bm_sum()`](https://inlabru-org.github.io/inlabru/reference/bm_sum.md),
[`bm_taylor()`](https://inlabru-org.github.io/inlabru/reference/bm_taylor.md),
[`bru_get_mapper()`](https://inlabru-org.github.io/inlabru/reference/bru_get_mapper.md),
[`bru_mapper()`](https://inlabru-org.github.io/inlabru/reference/bru_mapper.md)

## Examples

``` r
(m <- bm_collect(list(
  a = bm_index(2),
  b = bm_index(3)
), hidden = FALSE))
#> collect(a = index, b = index)
ibm_eval2(m, list(a = c(1, 2), b = c(1, 3, 2)), 1:5)
#> $offset
#> [1] 1 2 3 5 4
#> 
#> $jacobian
#> 5 x 5 sparse Matrix of class "dgTMatrix"
#>               
#> [1,] 1 . . . .
#> [2,] . 1 . . .
#> [3,] . . 1 . .
#> [4,] . . . . 1
#> [5,] . . . 1 .
#> 
```
