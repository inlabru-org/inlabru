# Mapper for linear Taylor approximations

Provides a pre-computed affine mapping, internally used to represent and
evaluate linearisation information. The `state0` information indicates
for which state the `offset` was evaluated; The affine mapper output is
defined as `effect(state) = offset + jacobian %*% (state - state0)`

## Usage

``` r
bm_taylor(offset = NULL, jacobian = NULL, state0 = NULL, values_mapper = NULL)

bru_mapper_taylor(...)
```

## Arguments

- offset:

  For `bm_taylor`, an offset vector giving the value of the
  linearisation at `state0`. May be `NULL`, interpreted as an all-zero
  vector of length determined by a non-null Jacobian.

- jacobian:

  For `bm_taylor()`, the Jacobian matrix, evaluated at `state0`, or, a
  named list of such matrices. May be `NULL` or an empty list, for a
  constant mapping.

- state0:

  For `bm_taylor`, the reference `state` for the linearisation, or a
  list of such states matching the `jacobian` list. `NULL` is
  interpreted as 0.

- values_mapper:

  mapper object to be used for `ibm_n` and `ibm_values` for
  `inla_f=TRUE` (experimental, currently unused)

- ...:

  Arguments passed on to `bm_taylor()`

## Value

A `bm_taylor` mapper object.

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
[`bm_repeat()`](https://inlabru-org.github.io/inlabru/reference/bm_repeat.md),
[`bm_scale()`](https://inlabru-org.github.io/inlabru/reference/bm_scale.md),
[`bm_shift()`](https://inlabru-org.github.io/inlabru/reference/bm_shift.md),
[`bm_sum()`](https://inlabru-org.github.io/inlabru/reference/bm_sum.md),
[`bru_get_mapper()`](https://inlabru-org.github.io/inlabru/reference/bru_get_mapper.md),
[`bru_mapper()`](https://inlabru-org.github.io/inlabru/reference/bru_mapper.md)

## Examples

``` r
m <- bm_taylor(
  offset = rep(2, 3),
  jacobian = matrix(1:6, 3, 2),
  state0 = c(1, 2)
)
ibm_eval2(m, state = 2:3)
#> $offset
#> [1]  7  9 11
#> 
#> $jacobian
#>      [,1] [,2]
#> [1,]    1    4
#> [2,]    2    5
#> [3,]    3    6
#> 
```
