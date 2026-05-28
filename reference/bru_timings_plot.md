# Plot inlabru iteration timings

Draws the time per iteration for preprocessing (including
linearisation), [`inla()`](https://rdrr.io/pkg/INLA/man/inla.html)
calls, and line search. Iteration `0` is the time used for defining the
model structure.

## Usage

``` r
bru_timings_plot(x)
```

## Arguments

- x:

  a [bru](https://inlabru-org.github.io/inlabru/reference/bru.md)
  object, typically a result from
  [`bru()`](https://inlabru-org.github.io/inlabru/reference/bru.md) for
  a nonlinear predictor model

## Details

Requires the "ggplot2" package to be installed.

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- bru(...)
bru_timings_plot(fit)
} # }
```
