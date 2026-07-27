# Mapper for aggregation

Constructs a mapper that aggregates elements of the input state, so it
can be used e.g. for weighted summation or integration over blocks of
values.

## Usage

``` r
bm_aggregate(rescale = FALSE, n_block = NULL, type = NULL)

bru_mapper_aggregate(...)
```

## Arguments

- rescale:

  logical; For `bm_aggregate` and `bm_logsumexp`, specifies if the
  blockwise sums should be normalised by the blockwise weight sums or
  not:

  - `FALSE`: (default) Straight weighted sum, no rescaling.

  - `TRUE`: Divide by the sum of the weight values within each block.
    This is useful for integration averages, when the given weights are
    plain integration weights. If the weights are `NULL` or all ones,
    this is the same as dividing by the number of entries in each block.

- n_block:

  Predetermined number of output blocks. If `NULL`, overrides the
  maximum block index in the inputs. The priority order is
  `input$n_block`, the mapper definition `n_block`, then
  `max(input$block)`.

- type:

  character; if non-NULL, overrides the `rescale` argument, and
  constructs an aggregation mapper of the given type instead. Supported
  values are "sum", "average" (for regular `bm_aggregate()`),
  "logsumexp", "logaverageexp" (for
  [`bm_logsumexp()`](https://inlabru-org.github.io/inlabru/reference/bm_logsumexp.md)),
  and "logitaverage" (for
  [`bm_logitaverage()`](https://inlabru-org.github.io/inlabru/reference/bm_logitaverage.md)).

- ...:

  Arguments passed on to `bm_aggregate()`

## Value

A `bm_aggregate` mapper object.

## See also

[bru_mapper](https://inlabru-org.github.io/inlabru/reference/bru_mapper.md),
[bru_mapper_generics](https://inlabru-org.github.io/inlabru/reference/bru_mapper_generics.md)

Other mappers:
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
[`bm_sum()`](https://inlabru-org.github.io/inlabru/reference/bm_sum.md),
[`bm_taylor()`](https://inlabru-org.github.io/inlabru/reference/bm_taylor.md),
[`bru_get_mapper()`](https://inlabru-org.github.io/inlabru/reference/bru_get_mapper.md),
[`bru_mapper()`](https://inlabru-org.github.io/inlabru/reference/bru_mapper.md)

## Examples

``` r
m <- bm_aggregate()
ibm_eval2(m, list(block = c(1, 2, 1, 2), weights = 1:4), 11:14)
#> $offset
#> [1] 50 80
#> 
#> $jacobian
#> 2 x 4 sparse Matrix of class "dgCMatrix"
#>             
#> [1,] 1 . 3 .
#> [2,] . 2 . 4
#> 
ibm_eval2(m, list(block = c(1, 2, 1, 2), weights = 1:4, n_block = 3), 11:14)
#> $offset
#> [1] 50 80  0
#> 
#> $jacobian
#> 3 x 4 sparse Matrix of class "dgCMatrix"
#>             
#> [1,] 1 . 3 .
#> [2,] . 2 . 4
#> [3,] . . . .
#> 
```
