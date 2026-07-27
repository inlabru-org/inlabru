# Extract block-averaged GCPO scores from one or more fitted bru models

After fitting a model with `control.gcpo = list(enable = TRUE)`, this
function reads the raw per-observation GCPO vector from `fit$gcpo$gcpo`
and returns one score per block or observation group.

Two cases are handled automatically:

- Block-based (lgcp / `"cp"` family):

  When `BRU_block` is present in `response_data` (populated
  automatically by the `"cp"` family machinery from the samplers
  structure), INLA repeats the same GCPO value for every observation in
  a leave-out block. The mean within each block collapses the repeated
  values to one score per block.

- Groups-based (all other families):

  When `BRU_block` is absent, INLA builds groups internally (via
  `num.level.sets` or `friends`) or from a user-supplied `groups` list.
  One score per observation is returned directly from `fit$gcpo$gcpo`.

Scores are returned on the probability scale, consistent with
`fit$cpo$cpo`. For multi-likelihood models, cumulative row offsets into
the stacked INLA response vector are handled automatically.

## Usage

``` r
bru_block_gcpo(fit, ...)
```

## Arguments

- fit:

  A fitted object of class `bru`, or a named list of such objects, or
  several such objects passed via `\dots`. Each must have been fitted
  with `control.gcpo = list(enable = TRUE)` so that `fit$gcpo$gcpo` is
  non-NULL.

- ...:

  Additional named fitted `bru` objects.

## Value

For a single fit, a list with two elements:

- `blocks`:

  Named list with one element per likelihood, each a vector of block
  labels (for `"cp"` models: unique `BRU_block` values; for other
  models: observation indices).

- `gcpo`:

  GCPO scores on the probability scale, consistent with `fit$cpo$cpo`. A
  numeric vector for single-likelihood models; a named list of numeric
  vectors for multi-likelihood models.

For multiple fits, a named list of such objects, one per fit.

## See also

[`bru()`](https://inlabru-org.github.io/inlabru/reference/bru.md),
[`bru_obs()`](https://inlabru-org.github.io/inlabru/reference/bru_obs.md),
[`bru_obs_control_gcpo()`](https://inlabru-org.github.io/inlabru/reference/bru_obs_gcpo.md),
[`bru_gcpo_table()`](https://inlabru-org.github.io/inlabru/reference/bru_gcpo_table.md)

## Examples

``` r
# \donttest{
if (bru_safe_inla()) {
  # Block-based example (lgcp)
  cvpart <- cv_hex(
    gorillas_sf$boundary,
    cellsize = 5, # Set to 0.5 for a more realistic, but slower, example
    n_group = 3,
    resolution = c(95, 80)
  )
  cvpart$block_ID <- seq_len(nrow(cvpart))
  cvpart$group <- NULL

  # Coarser mesh to make the example run faster:
  mesh <- fmesher::fm_mesh_2d(
    boundary = list(
      gorillas_sf$boundary,
      fmesher::fm_segm(gorillas_sf$mesh, boundary = TRUE)
    ),
    crs = fmesher::fm_crs(gorillas_sf$mesh)
  )

  # Use all nests for a more realistic, but slower, example
  nests <- gorillas_sf$nests[1:10, , drop = FALSE]
  a <- sf::st_intersects(nests, cvpart)
  nests$.block <- unlist(a)

  fit1 <- lgcp(
    geometry ~ Intercept(1),
    data = nests,
    samplers = cvpart,
    domain = list(geometry = mesh),
    control.gcpo = list(enable = TRUE, type.cv = "joint")
  )
  result <- bru_block_gcpo(fit1)
  result$blocks
  result$gcpo
  sum(log(result$gcpo))

  # Multiple models
  pcmatern <- INLA::inla.spde2.pcmatern(
    mesh,
    prior.sigma = c(1, 0.01),
    prior.range = c(0.1, 0.01)
  )
  fit2 <- lgcp(
    geometry ~ Intercept(1) + field(geometry, model = pcmatern),
    data = nests,
    samplers = cvpart,
    domain = list(geometry = mesh),
    control.gcpo = list(enable = TRUE, type.cv = "joint")
  )
  result <- bru_block_gcpo(list(Model1 = fit1, Model2 = fit2))
  str(result)
}
#> List of 2
#>  $ Model1:List of 2
#>   ..$ blocks:List of 1
#>   .. ..$ lhood1: int [1:4] 2 1 3 4
#>   ..$ gcpo  : num [1:4] 3.45e-07 1.11e-01 4.53e-02 7.06e-01
#>  $ Model2:List of 2
#>   ..$ blocks:List of 1
#>   .. ..$ lhood1: int [1:4] 2 1 3 4
#>   ..$ gcpo  : num [1:4] 6.56e-07 1.23e-01 5.79e-02 7.24e-01
# }
```
