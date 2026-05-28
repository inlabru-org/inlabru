# Plot method for posterior marginals estimated by bru

From version `2.11.0`, `plot.bru(x, ...)` calls `plot.inla(x, ...)` from
the `INLA` package, unless the first argument after `x` is a
`character`, in which case the pre-`2.11.0` behaviour is used, calling
`plotmarginal.inla(x, ...)` instead.

Requires the `ggplot2` package.

## Usage

``` r
# S3 method for class 'bru'
plot(x, ...)

plotmarginal.inla(
  result,
  varname = NULL,
  index = NULL,
  link = function(x) {
     x
 },
  add = FALSE,
  ggp = TRUE,
  lwd = 3,
  ...
)
```

## Arguments

- x:

  a fitted
  [`bru()`](https://inlabru-org.github.io/inlabru/reference/bru.md)
  model.

- ...:

  Options passed on to other methods.

- result:

  an `inla` or `bru` result object

- varname:

  character; name of the variable to plot

- index:

  integer; index of the random effect to plot

- link:

  function; link function to apply to the variable

- add:

  logical; if `TRUE`, add to an existing plot

- ggp:

  logical; unused

- lwd:

  numeric; line width

## Examples

``` r
if (FALSE) { # \dontrun{
if (require("ggplot2", quietly = TRUE)) {
  # Generate some data and fit a simple model
  input.df <- data.frame(x = cos(1:10)) |>
    dplyr::mutate(
      y = 5 + 2 * x + rnorm(length(x), mean = 0, sd = 0.1)
    )
  fit <- bru(y ~ x, family = "gaussian", data = input.df)
  summary(fit)

  # Plot the posterior density of the model's x-effect
  plot(fit, "x")
}
} # }
```
