# Mapper for cos/sin functions

Constructs a mapper for `cos`/`sin` functions of orders 1 (if
`intercept` is `TRUE`, otherwise 0) through `order`. The total number of
basis functions is `intercept + 2 * order`.

Optionally, each order can be given a non-unit scaling, via the
`scaling` vector, of length `intercept + order`. This can be used to
give an effective spectral prior. For example, let

    scaling = 1 / (1 + (0:4)^2)
    x <- seq(0, 1, length.out = 11)
    bmh1 = bm_harmonics(order = 4, interval = c(0, 1))
    u1 <- ibm_eval(
      bmh1,
      input = x,
      state = rnorm(9, sd = rep(scaling, c(1, 2, 2, 2, 2)))
    )

Then, with

    bmh2 = bm_harmonics(order = 4, scaling = scaling)
    u2 = ibm_eval(bmh2, input = x, state = rnorm(9))

the stochastic properties of `u1` and `u2` will be the same, with
`scaling^2` determining the variance for each frequency contribution.

The period for the first order harmonics is shifted and scaled to match
`interval`.

## Usage

``` r
bm_harmonics(order = 1, scaling = 1, intercept = TRUE, interval = c(0, 1))

bru_mapper_harmonics(...)
```

## Arguments

- order:

  For `bm_harmonics`, specifies the maximum `cos`/`sin` order. (Default
  1)

- scaling:

  For `bm_harmonics`, specifies an optional vector of scaling factors of
  length `intercept + order`, or a common single scalar.

- intercept:

  logical; For `bm_harmonics`, if `TRUE`, the first basis function is a
  constant. (Default `TRUE`)

- interval:

  numeric length-2 vector specifying a domain interval. Default
  `c(0, 1)`.

- ...:

  Arguments passed on to `bm_harmonics()`

## Value

A `bm_harmonics` mapper object.

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
m <- bm_harmonics(2)
ibm_eval2(m, input = c(0, pi / 4, pi / 2, 3 * pi / 4), 1:5)
#> $offset
#> [1]  7.000000 -7.247183  4.306719 -3.678570
#> 
#> $jacobian
#> 4 x 5 Matrix of class "dgeMatrix"
#>      [,1]       [,2]       [,3]       [,4]       [,5]
#> [1,]    1  1.0000000  0.0000000  1.0000000  0.0000000
#> [2,]    1  0.2205840 -0.9753680 -0.9026854 -0.4303012
#> [3,]    1 -0.9026854 -0.4303012  0.6296817  0.7768532
#> [4,]    1 -0.6188200  0.7855328 -0.2341236 -0.9722068
#> 
```
