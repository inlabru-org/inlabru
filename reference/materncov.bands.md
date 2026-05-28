# Matern correlation or covariance function approximate credible bands.

Evaluate the covariance function for an inla.spde object. Calculates the
posterior distribution of the range, log(range), variance, or
log(variance) parameter of a model's SPDE component. Can also calculate
Matern correlation or covariance function.

## Usage

``` r
materncov.bands(
  manifold,
  dist,
  log.range,
  log.variance = NULL,
  alpha = 2,
  quantile = 0.95,
  n = 64,
  S1.L = NULL
)
```

## Arguments

- manifold:

  Either "R1", "S1", "R2", or "S2", from `fm_manifold(mesh)`, or an
  object supported by
  [`fmesher::fm_manifold()`](https://inlabru-org.github.io/fmesher/reference/fm_manifold.html).

- dist:

  A vector of distances at which to calculate the
  covariances/correlations

- log.range:

  A scalar or a list (mean, sd), such as produced by
  `inla.spde.result(...)$summary.log.range.nominal[[1]][c("mean","sd")]`

- log.variance:

  Either `NULL`, a scalar, or vector of the same type as for log.range.
  When `NULL`, the correlations are calculated instead of the
  covariances.

- alpha:

  The SPDE operator order. Default 2.

- quantile:

  The target credible probability. Default 0.95.

- n:

  The number of parameter combinations to use for the approximation.
  Default 64.

- S1.L:

  For `manifold` `"S1"`, give the length of the cyclic interval

## Value

A list with estimated covariance or correlation (when `log.variance` is
`NULL`) functions:

- lower:

  An approximate lower bound for the `quantile` credible region

- median:

  The function for for the approximate median parameters quantile

- upper:

  An approximate upper bound for the `quantile` credible region

## Details

Uses a Gaussian assumption for the internal model parameters, and finds
a region in parameter space with approximately `quantile` probability.

## Author

Finn Lindgren <Finn.Lindgren@ed.ac.uk>
