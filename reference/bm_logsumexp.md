# Mapper for log-sum-exp aggregation

Constructs a mapper that aggregates elements of `exp(state)`, with
optional non-negative weighting, and then takes the
[`log()`](https://rdrr.io/r/base/Log.html), so it can be used e.g. for
\\v_k=\log\[\sum\_{i\in I_k} w_i \exp(u_i)\] \\ and
\\v_k=\log\[\sum\_{i\in I_k} w_i \exp(u_i) / \sum\_{i\in I_k} w_i\] \\
calculations. Relies on the input handling methods for `bm_aggregate`,
but also allows the weights to be supplied on a logarithmic scale as
`log_weights`. To avoid numerical overflow, it uses the common method of
internally shifting the state blockwise; \\v_k=s_k+\log\[\sum\_{i\in
I_k} \exp(u_i + \log(w_i)- s_k)\] \\, where \\s_k=\max\_{i\in I_k} u_i +
\log(w_i)\\ is the shift for block \\k\\.

## Usage

``` r
bm_logsumexp(rescale = FALSE, n_block = NULL)

bru_mapper_logsumexp(...)
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

- ...:

  Arguments passed on to `bm_logsumexp()`

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
m <- bm_logsumexp()
ibm_eval2(m, list(block = c(1, 2, 1, 2), weights = 1:4), 11:14)
#> $offset
#> [1] 14.14274 15.45177
#> 
#> $jacobian
#> 2 x 4 sparse Matrix of class "dgCMatrix"
#>                                               
#> [1,] 0.04316453 .          0.9568355 .        
#> [2,] .          0.06337894 .         0.9366211
#> 
ibm_eval2(m, list(block = c(1, 2, 1, 2), weights = 1:4, n_block = 3), 11:14)
#> $offset
#> [1] 14.14274 15.45177     -Inf
#> 
#> $jacobian
#> 3 x 4 sparse Matrix of class "dgCMatrix"
#>                                               
#> [1,] 0.04316453 .          0.9568355 .        
#> [2,] .          0.06337894 .         0.9366211
#> [3,] .          .          .         .        
#> 
```
