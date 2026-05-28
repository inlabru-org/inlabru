# Mapper for marginal distribution transformation

Constructs a mapper that transforms the marginal distribution `state`
from \\\textrm{N}(0,1)\\ to the distribution of a given (continuous)
quantile function. The `...` arguments are used as parameter arguments
to `qfun`, `pfun`, `dfun`, and `dqfun`.

## Usage

``` r
bm_marginal(qfun, pfun = NULL, dfun = NULL, dqfun = NULL, ..., inverse = FALSE)

bru_mapper_marginal(...)
```

## Arguments

- qfun:

  A quantile function, supporting arguments `p`, `lower.tail`, and
  `log.p`, like [`stats::qnorm()`](https://rdrr.io/r/stats/Normal.html).

- pfun:

  A CDF, supporting arguments `q`, `lower.tail`, and `log.p`, like
  [`stats::pnorm()`](https://rdrr.io/r/stats/Normal.html). Only needed
  and used when `xor(mapper[["inverse"]], reverse)` is `TRUE` in a
  method call. Default `NULL`

- dfun:

  A pdf, supporting arguments `x` and `log`, like
  [`stats::dnorm()`](https://rdrr.io/r/stats/Normal.html). If `NULL`
  (default), uses finite differences on `qfun` or `pfun` instead.

- dqfun:

  A function evaluating the reciprocal of the derivative of `qfun`, i.e.
  the density at a given CDF value. If `NULL` (default), uses
  `dfun(qfun(...),...)` or finite differences on `qfun` or `pfun`
  instead. Must support the same arguments as `qfun`.

- ...:

  Arguments passed on to `bm_marginal()`

- inverse:

  logical; If `FALSE` (default), `bm_marginal()` defines a mapping from
  standard Normal to a specified distribution. If `TRUE`, it defines a
  mapping from the specified distribution to a standard Normal.

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
m <- bm_marginal(qexp, pexp, rate = 1 / 8)
(val <- ibm_eval(m, state = -5:5))
#>  [1] 2.293213e-06 2.533739e-04 1.080648e-02 1.841033e-01 1.382030e+00
#>  [6] 5.545177e+00 1.472817e+01 3.026547e+01 5.286181e+01 8.288081e+01
#> [11] 1.205200e+02
ibm_eval(m, state = val, reverse = TRUE)
#>  [1] -5 -4 -3 -2 -1  0  1  2  3  4  5
m <- bm_marginal(qexp, pexp, dexp, rate = 1 / 8)
ibm_eval2(m, state = -3:3)
#> $offset
#> [1]  0.01080648  0.18410327  1.38203023  5.54517744 14.72817316 30.26547467
#> [7] 52.86180977
#> 
#> $jacobian
#> 7 x 7 diagonal matrix of class "ddiMatrix"
#>            [,1]      [,2]   [,3]     [,4]     [,5]     [,6]     [,7]
#> [1,] 0.03550271         .      .        .        .        .        .
#> [2,]          . 0.4419829      .        .        .        .        .
#> [3,]          .         . 2.3008        .        .        .        .
#> [4,]          .         .      . 6.383076        .        .        .
#> [5,]          .         .      .        . 12.20108        .        .
#> [6,]          .         .      .        .        . 18.98572        .
#> [7,]          .         .      .        .        .        . 26.26479
#> 
```
