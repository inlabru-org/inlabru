# Mapper for linking several mappers in sequence

Create a pipe mapper, where `mappers` is a list of mappers, and the
evaluated output of each mapper is handed as the state to the next
mapper. The `input` format for the `ibm_eval` and `ibm_jacobian` methods
is a list of inputs, one for each mapper.

## Usage

``` r
bm_pipe(mappers)

bru_mapper_pipe(...)
```

## Arguments

- mappers:

  A list of `bru_mapper` objects

- ...:

  Arguments passed on to `bm_pipe()`

## See also

[bru_mapper](https://inlabru-org.github.io/inlabru/reference/bru_mapper.md),
[bru_mapper_generics](https://inlabru-org.github.io/inlabru/reference/bru_mapper_generics.md)

Other mappers:
[`bm_aggregate()`](https://inlabru-org.github.io/inlabru/reference/bm_aggregate.md),
[`bm_collect()`](https://inlabru-org.github.io/inlabru/reference/bm_collect.md),
[`bm_const()`](https://inlabru-org.github.io/inlabru/reference/bm_const.md),
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
m <- bm_pipe(list(
  scale = bm_scale(),
  shift = bm_shift()
))
ibm_eval2(m, input = list(scale = 2, shift = 1:4), state = 1:4)
#> $offset
#> [1]  3  6  9 12
#> 
#> $jacobian
#> 4 x 4 diagonal matrix of class "ddiMatrix"
#>      [,1] [,2] [,3] [,4]
#> [1,]    2    .    .    .
#> [2,]    .    2    .    .
#> [3,]    .    .    2    .
#> [4,]    .    .    .    2
#> 
```
