# Convenient model fitting using (iterated) INLA

This method is a wrapper for
[`INLA::inla`](https://rdrr.io/pkg/INLA/man/inla.html) and provides
multiple enhancements.

- Easy usage of spatial covariates and automatic construction of inla
  projection matrices for (spatial) SPDE models. This feature is
  accessible via the `components` parameter. Practical examples on how
  to use spatial data by means of the components parameter can also be
  found by looking at the
  [`lgcp()`](https://inlabru-org.github.io/inlabru/reference/lgcp.md)
  function's documentation.

- Constructing multiple observation models is straightforward. See
  [`bru_obs()`](https://inlabru-org.github.io/inlabru/reference/bru_obs.md)
  for more information on how to provide additional models to `bru`
  using the `...` parameter list.

- Support for non-linear predictors. See example below.

- Log Gaussian Cox process (LGCP) inference is available by using the
  `"cp"` family or (even easier) by using the
  [`lgcp()`](https://inlabru-org.github.io/inlabru/reference/lgcp.md)
  function.

## Usage

``` r
bru(components = ~Intercept(1), ..., options = list(), .envir = parent.frame())

bru_rerun(result, options = list())

# S3 method for class 'bru'
summary(object, verbose = FALSE, ...)

# S3 method for class 'summary_bru'
print(x, ...)

# S3 method for class 'bru'
print(x, ...)
```

## Arguments

- components:

  Latent component definitions, either as a
  [`bru_comp_list()`](https://inlabru-org.github.io/inlabru/reference/bru_comp_list.md)
  object, or a `formula`-like specification. Also used to define a
  default linear additive predictor. See
  [`bru_comp()`](https://inlabru-org.github.io/inlabru/reference/bru_comp.md)
  for details.

- ...:

  Observation models, each constructed by a calling
  [`bru_obs()`](https://inlabru-org.github.io/inlabru/reference/bru_obs.md),
  or
  [`bru_obs_list()`](https://inlabru-org.github.io/inlabru/reference/bru_obs.md).

  Alternatively, for backwards compatibility, may be named parameters
  that can be passed to a single
  [`bru_obs()`](https://inlabru-org.github.io/inlabru/reference/bru_obs.md)
  call. These arguments will be evaluated before calling
  [`bru_obs()`](https://inlabru-org.github.io/inlabru/reference/bru_obs.md),
  in order to detect if they already are `bru_obs` objects. This means
  that special arguments that are only available in the context of
  `data` or `response_data` (such as `Ntrials`) will only work properly
  in direct calls to
  [`bru_obs()`](https://inlabru-org.github.io/inlabru/reference/bru_obs.md).

- options:

  A
  [bru_options](https://inlabru-org.github.io/inlabru/reference/bru_options.md)
  options object or a list of options passed on to
  [`bru_options()`](https://inlabru-org.github.io/inlabru/reference/bru_options.md)

- .envir:

  Environment for component evaluation (for when a non-formula
  specification is used)

- result:

  A previous estimation object of class `bru`

- object:

  An object obtained from a `bru()` or
  [`lgcp()`](https://inlabru-org.github.io/inlabru/reference/lgcp.md)
  call

- verbose:

  logical; If `TRUE`, include more details of the component definitions.
  If `FALSE`, only show basic component definition information. Default:
  `FALSE`

- x:

  An object to be printed

## Value

bru returns an object of class "bru". A `bru` object inherits from
[`INLA::inla`](https://rdrr.io/pkg/INLA/man/inla.html) (see the inla
documentation for its properties) and adds additional information stored
in the `bru_info` field.

## Methods (by generic)

- `summary(bru)`: Takes a fitted `bru` object produced by `bru()` or
  [`lgcp()`](https://inlabru-org.github.io/inlabru/reference/lgcp.md)
  and creates various summaries from it, including the summary output
  from the corresponding `INLA::sumary()` method. The ... arguments are
  passed on to component summary functions, see
  [`summary.bru_comp()`](https://inlabru-org.github.io/inlabru/reference/summary.bru_comp.md).

- `print(bru)`: Print a summary of a `bru` object.

## Functions

- `bru_rerun()`: Continue the optimisation from a previously computed
  estimate. The estimation `options` list can be given new values to
  override the original settings.

  To rerun with a subset of the data (e.g. for cross validation or prior
  sampling), use
  [`bru_set_missing()`](https://inlabru-org.github.io/inlabru/reference/bru_set_missing.md)
  to set all or part of the response data to `NA` before calling
  `bru_rerun()`.

## Author

Fabian E. Bachl <bachlfab@gmail.com>

## Examples

``` r
# \donttest{
if (bru_safe_inla()) {
  # Simulate some covariates x and observations y
  input.df <- data.frame(x = cos(1:10))
  input.df <- within(input.df, {
    y <- 5 + 2 * x + rnorm(10, mean = 0, sd = 0.1)
  })

  # Fit a Gaussian likelihood model
  fit <- bru(y ~ x + Intercept(1), family = "gaussian", data = input.df)

  # Obtain summary
  fit$summary.fixed
}
#> Changing INLA option num.threads from '4:1' to '1:1:1'.
#>               mean        sd 0.025quant 0.5quant 0.975quant     mode
#> x         2.027702 0.0521093   1.923574 2.027703   2.131825 2.027703
#> Intercept 4.959784 0.0368416   4.886164 4.959785   5.033399 4.959784
#>                    kld
#> x         5.770000e-06
#> Intercept 5.770205e-06


if (bru_safe_inla()) {
  # Alternatively, we can use the bru_obs() function to construct the
  # likelihood:

  lik <- bru_obs(
    family = "gaussian",
    formula = y ~ x + Intercept,
    data = input.df
  )
  fit <- bru(~ x + Intercept(1), lik)
  fit$summary.fixed
}
#>               mean        sd 0.025quant 0.5quant 0.975quant     mode
#> x         2.027702 0.0521093   1.923574 2.027703   2.131825 2.027703
#> Intercept 4.959784 0.0368416   4.886164 4.959785   5.033399 4.959784
#>                    kld
#> x         5.770000e-06
#> Intercept 5.770205e-06

# An important addition to the INLA methodology is bru's ability to use
# non-linear predictors. Such a predictor can be formulated via bru_obs()'s
# \code{formula} parameter. The z(1) notation is needed to ensure that
# the z component should be interpreted as single latent variable and not
# a covariate:

if (bru_safe_inla()) {
  z <- 2
  input.df <- within(input.df, {
    y <- 5 + exp(z) * x + rnorm(10, mean = 0, sd = 0.1)
  })
  lik <- bru_obs(
    family = "gaussian", data = input.df,
    formula = y ~ exp(z) * x + Intercept
  )
  fit <- bru(~ z(1) + Intercept(1), lik)

  # Check the result (z posterior should be around 2)
  fit$summary.fixed
}
#>               mean         sd 0.025quant 0.5quant 0.975quant     mode
#> z         1.995419 0.00791756   1.979598 1.995419    2.01124 1.995419
#> Intercept 5.021860 0.04117309   4.939585 5.021861    5.10413 5.021861
#>                    kld
#> z         5.753208e-06
#> Intercept 5.752969e-06
# }
```
