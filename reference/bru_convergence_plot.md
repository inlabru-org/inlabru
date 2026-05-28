# Plot inlabru convergence diagnostics

Draws four panels of convergence diagnostics for an iterated INLA method
estimation

## Usage

``` r
bru_convergence_plot(x, from = 1, to = NULL, type = NULL)
```

## Arguments

- x:

  a [bru](https://inlabru-org.github.io/inlabru/reference/bru.md)
  object, typically a result from
  [`bru()`](https://inlabru-org.github.io/inlabru/reference/bru.md) for
  a nonlinear predictor model

- from, to:

  integer values for the range of iterations to plot. Default `from = 1`
  (start from the first iteration) and `to = NULL` (end at the last
  iteration). Set `from = 0` to include the initial linearisation point
  in the track plot.

- type:

  **\[experimental\]** character; "bru" (default) for iterative
  nonlinear inlabru convergence diagnostics plots, or "inla" for INLA
  optimiser trace plots.

## Value

A ggplot object with four panels of convergence diagnostics:

- `Tracks`: Mode and linearisation values for each effect

- `Mode - Lin`: Difference between mode and linearisation values for
  each effect

- `|Change| / sd`: Absolute change in mode and linearisation values
  divided by the standard deviation for each effect

- `Change & sd`: Absolute change in mode and linearisation values and
  standard deviation for each effect

For multidimensional components, only the overall average, maximum, and
minimum values are shown.

## Details

Requires the "dplyr", "ggplot2", and "patchwork" packages to be
installed.

## See also

[`bru()`](https://inlabru-org.github.io/inlabru/reference/bru.md)

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- bru(...)
bru_convergence_plot(fit)
} # }
```
