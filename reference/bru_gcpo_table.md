# Compare block-averaged GCPO scores across multiple fitted bru models

Calls
[`bru_block_gcpo()`](https://inlabru-org.github.io/inlabru/reference/bru_block_gcpo.md)
on each fit and combines the results into a `data.frame` with one row
per block and one column per model, making it straightforward to compare
GCPO scores across models block by block.

## Usage

``` r
bru_gcpo_table(fits = NULL, ...)
```

## Arguments

- fits:

  A named list of fitted `bru` objects, or `NULL` if models are passed
  via `\dots`.

- ...:

  Named fitted `bru` objects, used when `fits` is not a list.

## Value

For a single-likelihood model, a `data.frame` with columns:

- `block`:

  Block labels from `response_data$BRU_block` (for `"cp"` models) or
  observation indices (for other models).

- one column per model:

  GCPO scores on the probability scale, named after the elements of
  `fits` or `\dots`.

For multi-likelihood models, a named list of such `data.frame`s, one per
likelihood.

## See also

[`bru_block_gcpo()`](https://inlabru-org.github.io/inlabru/reference/bru_block_gcpo.md)

## Examples

``` r
# \donttest{
if (bru_safe_inla()) {
  cvpart <- cv_hex(
    gorillas_sf$boundary,
    cellsize = 0.5,
    n_group = 3,
    resolution = c(95, 80)
  )
  cvpart$block_ID <- seq_len(nrow(cvpart))
  cvpart$group <- NULL

  nests <- gorillas_sf$nests
  a <- sf::st_intersects(nests, cvpart)
  nests$.block <- unlist(a)

  fit1 <- lgcp(
    geometry ~ Intercept(1),
    data = nests,
    samplers = cvpart,
    domain = list(geometry = gorillas_sf$mesh),
    control.gcpo = list(enable = TRUE, type.cv = "joint")
  )
  pcmatern <- INLA::inla.spde2.pcmatern(
    gorillas_sf$mesh,
    prior.sigma = c(1, 0.01),
    prior.range = c(0.1, 0.01)
  )
  fit2 <- lgcp(
    geometry ~ Intercept(1) + field(geometry, model = pcmatern),
    data = nests,
    samplers = cvpart,
    domain = list(geometry = gorillas_sf$mesh),
    control.gcpo = list(enable = TRUE, type.cv = "joint")
  )
  gcpo_df <- bru_gcpo_table(Model1 = fit1, Model2 = fit2)
  names(gcpo_df)
  colSums(log(gcpo_df[, -1]))
}
#>   Model1   Model2 
#> 1593.967 2395.980 
# }
```
