# Mapper for factor variables

Create a factor mapper

## Usage

``` r
bm_factor(values, factor_mapping, indexed = FALSE)

bru_mapper_factor(...)
```

## Arguments

- values:

  Input values calculated by
  [`bru_input.bru_input()`](https://inlabru-org.github.io/inlabru/reference/bru_input.md)

- factor_mapping:

  character; selects the type of factor mapping.

  - `'contrast'` for leaving out the first factor level.

  - `'full'` for keeping all levels.

- indexed:

  logical; if `TRUE`, the
  [`ibm_values()`](https://inlabru-org.github.io/inlabru/reference/ibm_values.md)
  method will return an integer vector instead of the factor levels.
  This is needed e.g. for `group` and `replicate` mappers, since
  [`INLA::f()`](https://rdrr.io/pkg/INLA/man/f.html) doesn't accept
  factor values. Default: `FALSE`, which works for the main input
  mappers. The default mapper constructions will set it the required
  setting.

- ...:

  Arguments passed on to `bm_factor()`

## Value

A `bm_factor` mapper object.

## See also

[bru_mapper](https://inlabru-org.github.io/inlabru/reference/bru_mapper.md),
[bru_mapper_generics](https://inlabru-org.github.io/inlabru/reference/bru_mapper_generics.md)

Other mappers:
[`bm_aggregate()`](https://inlabru-org.github.io/inlabru/reference/bm_aggregate.md),
[`bm_collect()`](https://inlabru-org.github.io/inlabru/reference/bm_collect.md),
[`bm_const()`](https://inlabru-org.github.io/inlabru/reference/bm_const.md),
[`bm_expr()`](https://inlabru-org.github.io/inlabru/reference/bm_expr.md),
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
m <- bm_factor(factor(c("a", "b")), "full")
ibm_eval2(m, input = c("b", "a", "a", "b"), state = c(1, 3))
#> $offset
#> [1] 3 1 1 3
#> 
#> $jacobian
#> 4 x 2 sparse Matrix of class "dgCMatrix"
#>         
#> [1,] . 1
#> [2,] 1 .
#> [3,] 1 .
#> [4,] . 1
#> 

m <- bm_factor(factor(c("a", "b")), "contrast")
ibm_eval2(m, input = factor(c("b", "a", "a", "b")), state = 2)
#> $offset
#> [1] 2 0 0 2
#> 
#> $jacobian
#> 4 x 1 sparse Matrix of class "dgCMatrix"
#>       
#> [1,] 1
#> [2,] .
#> [3,] .
#> [4,] 1
#> 
```
