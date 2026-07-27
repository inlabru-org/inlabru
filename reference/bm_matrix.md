# Mapper for matrix multiplication

Create a matrix mapper, for a given number of columns

## Usage

``` r
bm_matrix(labels)

bru_mapper_matrix(...)
```

## Arguments

- labels:

  Column labels for matrix mappings; Can be factor, character, or a
  single integer specifying the number of columns for integer column
  indexing.

- ...:

  Arguments passed on to `bm_matrix()`

## Value

A `bm_matrix` mapper object.

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
m <- bm_matrix(labels = c("a", "b"))
ibm_values(m)
#> [1] a b
#> Levels: a b
ibm_eval2(m, input = matrix(1:6, 3, 2), state = 2:3)
#> $offset
#> [1] 14 19 24
#> 
#> $jacobian
#> 3 x 2 Matrix of class "dgeMatrix"
#>      a b
#> [1,] 1 4
#> [2,] 2 5
#> [3,] 3 6
#> 

m <- bm_matrix(labels = 2L)
ibm_values(m)
#> [1] 1 2
ibm_eval2(m, input = matrix(1:6, 3, 2), state = 2:3)
#> $offset
#> [1] 14 19 24
#> 
#> $jacobian
#> 3 x 2 Matrix of class "dgeMatrix"
#>      1 2
#> [1,] 1 4
#> [2,] 2 5
#> [3,] 3 6
#> 
```
