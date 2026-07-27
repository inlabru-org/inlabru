# Mapper for `fm_mesh_1d`

Create mapper for an `fm_mesh_1d` object

## Usage

``` r
# S3 method for class 'fm_mesh_1d'
bru_mapper(mesh, indexed = TRUE, ...)
```

## Arguments

- mesh:

  An `fm_mesh_1d` object to use as a mapper

- indexed:

  logical; If `TRUE` (default), the
  [`ibm_values()`](https://inlabru-org.github.io/inlabru/reference/ibm_values.md)
  output will be the integer indexing sequence for the latent variables
  (needed for `spde` models). If `FALSE`, points representative of the
  basis centres are returned (useful for an interpolator for `rw2`
  models and similar, for `fmesher` versions `>= 0.3.0.9002`).

- ...:

  Arguments passed on to
  [`bm_fmesher()`](https://inlabru-org.github.io/inlabru/reference/bm_fmesher.md)

## Value

A `bm_fm_mesh_1d` or `bm_fmesher` object. The the general
[`bm_fmesher()`](https://inlabru-org.github.io/inlabru/reference/bm_fmesher.md)
mapper handles all indexed `fmesher` objects.

## See also

[bru_mapper](https://inlabru-org.github.io/inlabru/reference/bru_mapper.md),
[bru_mapper_generics](https://inlabru-org.github.io/inlabru/reference/bru_mapper_generics.md)

Other mappers:
[`bm_aggregate()`](https://inlabru-org.github.io/inlabru/reference/bm_aggregate.md),
[`bm_collect()`](https://inlabru-org.github.io/inlabru/reference/bm_collect.md),
[`bm_const()`](https://inlabru-org.github.io/inlabru/reference/bm_const.md),
[`bm_expr()`](https://inlabru-org.github.io/inlabru/reference/bm_expr.md),
[`bm_factor()`](https://inlabru-org.github.io/inlabru/reference/bm_factor.md),
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
m <- bru_mapper(fmesher::fm_mesh_1d(c(1:3, 5, 7)))
ibm_values(m)
#> [1] 1 2 3 4 5
ibm_eval(m, 1:7, 1:5)
#> [1] 1.0 2.0 3.0 3.5 4.0 4.5 5.0

m <- bru_mapper(fmesher::fm_mesh_1d(c(1:3, 5, 7)), indexed = FALSE)
ibm_values(m)
#> [1] 1 2 3 5 7
ibm_eval(m, 1:7, 1:5)
#> [1] 1.0 2.0 3.0 3.5 4.0 4.5 5.0

m <- bru_mapper(
  fmesher::fm_mesh_1d(c(1:3, 5, 7), degree = 2, boundary = "free"),
  indexed = FALSE
)
ibm_values(m)
#> [1] 0.5 1.5 2.5 4.0 6.0 8.0
ibm_eval(m, 1:7, 1:6)
#> [1] 1.500000 2.500000 3.333333 3.958333 4.500000 5.000000 5.500000
```
