# Mapper for general `fmesher` function space objects

Creates a mapper for general `fmesher` function space objects.

## Usage

``` r
bm_fmesher(mesh)

bru_mapper_fmesher(...)

# S3 method for class 'fm_mesh_2d'
bru_mapper(mesh, ...)
```

## Arguments

- mesh:

  An `fmesher` object to map, supported by
  [fmesher::fm_basis](https://inlabru-org.github.io/fmesher/reference/fm_basis.html)`(mesh, input)`
  and
  [fmesher::fm_dof](https://inlabru-org.github.io/fmesher/reference/fm_dof.html)`(mesh)`.

- ...:

  Arguments passed on to `bm_fmesher()`

## Value

A `bm_fmesher` object.

## Details

Handles indexed mapping for all `fmesher` classes that support
[`fm_dof()`](https://inlabru-org.github.io/fmesher/reference/fm_dof.html)
and
[`fm_basis()`](https://inlabru-org.github.io/fmesher/reference/fm_basis.html)
methods. For non-indexed mapping of `fm_mesh_1d` objects, use
`bru_mapper(mesh, indexed = FALSE)` which invokes the
[`bru_mapper.fm_mesh_1d()`](https://inlabru-org.github.io/inlabru/reference/bm_fm_mesh_1d.md)
method.

## Functions

- `bru_mapper(fm_mesh_2d)`: Equivalent to calling `bm_fmesher()`. Note:
  Prior to version `2.12.0.9021`, this returned a
  `bru_mapper_fm_mesh_2d` object. Also see the note for
  [`bru_mapper.fm_mesh_1d()`](https://inlabru-org.github.io/inlabru/reference/bm_fm_mesh_1d.md).

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
m <- bm_fmesher(fmesher::fmexample$mesh)
ibm_n(m)
#> [1] 292
ibm_eval(m, as.matrix(expand.grid(-2:2, -2:2)), seq_len(ibm_n(m)))
#>  [1] 266.64363 113.59173  56.61537 106.69945 256.01042 136.21372 150.83769
#>  [8] 152.68096 178.13683  59.93549 131.57928 163.75497 152.42684 182.22284
#> [15] 163.25104 149.28159 182.92309 137.14656 152.31003 276.18942 102.82566
#> [22] 134.25713 112.48311 112.19143 221.34490
```
