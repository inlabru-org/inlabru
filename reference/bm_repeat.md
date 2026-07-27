# Mapper for repeating a mapper

Defines a repeated-space mapper that sums the contributions for each
copy. The
[`ibm_n()`](https://inlabru-org.github.io/inlabru/reference/ibm_n.md)
method returns `ibm_n(mapper) * n_rep`, and
[`ibm_values()`](https://inlabru-org.github.io/inlabru/reference/ibm_values.md)
returns `seq_len(ibm_n(mapper))`.

## Usage

``` r
bm_repeat(mapper, n_rep, interleaved = FALSE)

bru_mapper_repeat(...)
```

## Arguments

- mapper:

  The mapper to be repeated.

- n_rep:

  The number of times to repeat the mapper. If a vector, the
  non-interleaved repeats are combined into a single repeat mapping, and
  combined with interleaved repeats via a
  [`bm_sum()`](https://inlabru-org.github.io/inlabru/reference/bm_sum.md)
  of mappers.

- interleaved:

  logical; if `TRUE`, the repeated mapping columns are interleaved;
  `(x1[1], x2[1], ..., x1[2], x2[2], ...)`. If `FALSE` (default), the
  repeated mapping columns are contiguous,
  `(x1[1], x1[2], ..., x2[1], x2[2], ...)`, and the Jacobian is a
  [`cbind()`](https://rdrr.io/r/base/cbind.html) of the Jacobians of the
  repeated mappers.

  If `n_rep` is a vector, `interleaved` should either be a single
  logical, or a vector of the same length. Each element applies to the
  corresponding `n_rep` repetition specification.

- ...:

  Arguments passed on to
  [`bm_scale()`](https://inlabru-org.github.io/inlabru/reference/bm_scale.md)

## Value

A `bm_repeat` or `bm_sum` object, or the original input `mapper`.

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
[`bm_scale()`](https://inlabru-org.github.io/inlabru/reference/bm_scale.md),
[`bm_shift()`](https://inlabru-org.github.io/inlabru/reference/bm_shift.md),
[`bm_sum()`](https://inlabru-org.github.io/inlabru/reference/bm_sum.md),
[`bm_taylor()`](https://inlabru-org.github.io/inlabru/reference/bm_taylor.md),
[`bru_get_mapper()`](https://inlabru-org.github.io/inlabru/reference/bru_get_mapper.md),
[`bru_mapper()`](https://inlabru-org.github.io/inlabru/reference/bru_mapper.md)

## Examples

``` r
(m0 <- bm_index(3))
#> index
(m <- bm_repeat(m0, 5))
#> repeat(5 x index)
ibm_n(m)
#> [1] 15
ibm_values(m)
#>  [1]  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15
ibm_jacobian(m, 1:3)
#> 3 x 15 sparse Matrix of class "dgCMatrix"
#>                                   
#> [1,] 1 . . 1 . . 1 . . 1 . . 1 . .
#> [2,] . 1 . . 1 . . 1 . . 1 . . 1 .
#> [3,] . . 1 . . 1 . . 1 . . 1 . . 1
ibm_eval(m, 1:3, seq_len(ibm_n(m)))
#> [1] 35 40 45

# Interleaving and grouping
(m <- bm_repeat(m0, c(2, 1, 2), c(TRUE, FALSE, FALSE)))
#> sum(1 = {repeat(2 x index, interleaved)}, 2 = {repeat(3 x index)})
ibm_n(m)
#> [1] 15
ibm_values(m)
#>  [1]  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15
ibm_jacobian(m, 1:3)
#> 3 x 15 sparse Matrix of class "dgCMatrix"
#>                                   
#> [1,] 1 1 . . . . 1 . . 1 . . 1 . .
#> [2,] . . 1 1 . . . 1 . . 1 . . 1 .
#> [3,] . . . . 1 1 . . 1 . . 1 . . 1
ibm_eval(m, 1:3, seq_len(ibm_n(m)))
#> [1] 33 40 47
```
