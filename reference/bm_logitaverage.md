# Mapper for logit-sum-inverse-logit aggregation

**\[experimental\]** Constructs a mapper that averages elements of
`plogis(state)`, with optional non-negative weighting, and then takes
the [`qlogis()`](https://rdrr.io/r/stats/Logistic.html). Relies on the
input handling methods for `bm_aggregate`. To avoid numerical issues, it
uses `plogis(x, log.p = TRUE)`, `plogis(-x, log.p = TRUE)`, and the
equivalent of two applications of
[`bm_logsumexp()`](https://inlabru-org.github.io/inlabru/reference/bm_logsumexp.md)
to evaluate \\\log(p_k)-\log(1-p_k)\\, where

\\p_k=\sum\_{i\in I_k} w_i / (1+e^{-\eta_i}) / \sum\_{i\in I_k} w_i \\

## Usage

``` r
bm_logitaverage(n_block = NULL)
```

## Arguments

- n_block:

  Predetermined number of output blocks. If `NULL`, overrides the
  maximum block index in the inputs. The priority order is
  `input$n_block`, the mapper definition `n_block`, then
  `max(input$block)`.

## Value

A `bm_logitaverage/bm_aggregate` mapper object.

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
m <- bm_logitaverage()
ibm_eval2(m, list(block = c(1, 2, 1, 2), weights = 1:4), 11:14)
#> $offset
#> [1] 12.04555 12.85907
#> 
#> $jacobian
#> 2 x 4 sparse Matrix of class "dgCMatrix"
#>                                             
#> [1,] 0.7112239 .         0.2887761 .        
#> [2,] .         0.7869824 .         0.2130176
#> 
ibm_eval2(m, list(block = c(1, 2, 1, 2), weights = 1:4, n_block = 3), 11:14)
#> $offset
#> [1] 12.04555 12.85907      NaN
#> 
#> $jacobian
#> 3 x 4 sparse Matrix of class "dgCMatrix"
#>                                             
#> [1,] 0.7112239 .         0.2887761 .        
#> [2,] .         0.7869824 .         0.2130176
#> [3,] .         .         .         .        
#> 
```
