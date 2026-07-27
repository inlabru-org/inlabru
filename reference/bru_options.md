# Create or update an options objects

Create a new options object, or merge information from several objects.

The `_get`, `_set`, and `_reset` functions operate on a global package
options override object. In many cases, setting options in specific
calls to
[`bru()`](https://inlabru-org.github.io/inlabru/reference/bru.md) is
recommended instead.

## Usage

``` r
bru_options(...)

as.bru_options(x = NULL)

bru_options_default()

bru_options_check(options, ignore_null = TRUE)

bru_options_get(name = NULL, include_default = TRUE)

bru_options_set(..., .reset = FALSE)

bru_options_reset()

bru_options_set_local(..., .reset = FALSE, .envir = parent.frame())
```

## Arguments

- ...:

  A collection of named options, optionally including one or more
  `bru_options` objects. Options specified later override the previous
  options.

- x:

  An object to be converted to an `bru_options` object.

- options:

  An `bru_options` object to be checked

- ignore_null:

  Ignore missing or NULL options.

- name:

  Either `NULL`, or single option name string, or character vector or
  list with option names, Default: NULL

- include_default:

  logical; If `TRUE`, the default options are included together with the
  global override options. Default: `TRUE`

- .reset:

  For `bru_options_set`, logical indicating if the global override
  options list should be emptied before setting the new option(s).

- .envir:

  The environment in which to set the options. Default:
  [`parent.frame()`](https://rdrr.io/r/base/sys.parent.html)

## Value

`bru_options()` returns a `bru_options` object.

For `as.bru_options()`, `NULL` or no input returns an empty
`bru_options` object, a `list` is converted via `bru_options(...)`, and
`bru_options` input is passed through. Other types of input generates an
error.

`bru_options_default()` returns an `bru_options` object containing
default options.

`bru_options_check()` returns a `logical`; `TRUE` if the object contains
valid options for use by other functions

`bru_options_get` returns either an `bru_options` object, for
`name == NULL`, the contents of single option, if `name` is a options
name string, or a named list of option contents, if `name` is a list of
option name strings.

`bru_options_set()` returns a copy of the global override options,
invisibly (as `bru_options_get(include_default = FALSE)`).

## Functions

- `as.bru_options()`: Coerces inputs to a `bru_options` object.

- `bru_options_default()`: Returns the default options.

- `bru_options_check()`: Checks for valid contents of a `bru_options`
  object, and produces warnings for invalid options.

- `bru_options_get()`: Used to access global package options.

- `bru_options_set()`: Used to set global package options.

- `bru_options_reset()`: Clears the global option overrides.

- `bru_options_set_local()`: Sets local option overrides, that are
  automatically reset using
  [`withr::defer()`](https://withr.r-lib.org/reference/defer.html).

## Valid options

For `bru_options` and `bru_options_set`, recognised options are:

- bru_verbose:

  numeric or logical; if `TRUE`, log messages with `verbosity` \\\le\\ 1
  are printed by
  [`bru_log_message()`](https://inlabru-org.github.io/inlabru/reference/bru_log_message.md).
  If numeric, messages with `verbosity` \\\le\\ `bru_verbose` are
  printed. For line search details, set `bru_verbose=2` or `3`. Default:
  0, to not print any messages.

- bru_verbose_store:

  numeric or logical, as for `bru_verbose`, but controls what messages
  are stored in the global log object. Default: `Inf`, to store all
  messages.

- bru_max_iter:

  maximum number of inla iterations, default 10. Also see
  `bru_method$rel_tol` and related options below.

- bru_run:

  If `TRUE`, run inference. Otherwise only return configuration needed
  to run inference.

- bru_initial:

  An `inla` object returned from previous calls of
  [`INLA::inla`](https://rdrr.io/pkg/INLA/man/inla.html),
  [`bru()`](https://inlabru-org.github.io/inlabru/reference/bru.md) or
  [`lgcp()`](https://inlabru-org.github.io/inlabru/reference/lgcp.md),
  or a list of named vectors of starting values for the latent
  variables. This will be used as a starting point for further
  improvement of the approximate posterior.

- bru_int_args:

  List of arguments passed all the way to the integration method
  [`fmesher::fm_int()`](https://inlabru-org.github.io/fmesher/reference/fm_int.html)
  for 'cp' family models;

  method

  :   "stable" or "direct". For "stable" (default) integration points
      are aggregated to mesh vertices.

  nsub1

  :   Number of integration points per knot interval in 1D. Default 30.

  nsub2

  :   Number of integration points along a triangle edge for 2D. Default
      9.

  nsub

  :   Deprecated parameter that overrides `nsub1` and `nsub2` if set.
      Default `<not set>`.

- bru_method:

  List of arguments controlling the iterative inlabru method:

  autodiff

  :   Controls the linearisation calculation method. One of

      'pandemic'

      :   Treats all post-component transformations as a single
          calculation (the default from version `2.1.15`).

      'fullchain'

      :   Uses the chain rule for the composition of components,
          predictor, and post-predictor transformations (new from
          `2.14.1.9008`).

  agg

  :   Controls the aggregation implementation. One of

      'pandemic'

      :   Treats all post-component transformations as a single
          calculation (the default from version `2.1.15`).

      'fullchain'

      :   Uses the chain rule for the composition of components,
          predictor, and post-predictor transformations (new from
          `2.14.1.9008`).

  finite_diff

  :   One of

      'forward'

      :   Use onesided finite differences.

      'central'

      :   Use symmetric differences.

  taylor

  :   'pandemic' (the default from version `2.1.15`).

  search

  :   Either 'all' (default), to use all available line search methods,
      or one or more of

      'finite'

      :   (reduce step size until predictor is finite)

      'contract'

      :   (decrease step size until trust hypersphere reached)

      'expand'

      :   (increase step size until no improvement)

      'optimise'

      :   (fast approximate error norm minimisation)

      To disable line search, set to an empty vector.

  factor

  :   Numeric, \\\> 1\\ determining the line search step scaling
      multiplier. Default \\(1 + \sqrt{5})/2\\.

  rel_tol

  :   Stop the iterations when the largest change in linearisation point
      (the conditional latent state mode) in relation to the estimated
      posterior standard deviation is less than `rel_tol`. Default 0.1
      (ten percent).

  max_step

  :   The largest allowed line search step factor. Factor 1 is the full
      INLA step. Default is 2.

  line_opt_method

  :   Which method to use for the line search optimisation step. Default
      "onestep", using a quadratic approximation based on the value and
      gradient at zero, and the value at the current best step length
      guess. The method "full" does line optimisation on the full
      nonlinear predictor; this is slow and intended for debugging
      purposes only.

- bru_compress_cp:

  logical; when `TRUE`, compress the \\\sum\_{i=1}^n \eta_i\\ part of
  the Poisson process likelihood (`family = "cp"`) into either a single
  term, with \\y=n\\, and predictor `mean(eta)`, or a blockwise version
  of this. Default: `FALSE` (was `TRUE` prior to version `2.13.0.9031`.)

- bru_debug:

  logical; when `TRUE`, activate temporary debug features for package
  development. Default: `FALSE`

- bru_compat_pre_2_14_enable:

  logical; when `TRUE`, enable compatibility features for inlabru
  versions prior to 2.14. Set to `FALSE` to test external package
  compatibility updates. Default: `TRUE` before version 2.14, and will
  be set to `FALSE` by default in a later version.

- [`inla()`](https://rdrr.io/pkg/INLA/man/inla.html) options:

  All options not starting with `bru_` are passed on to
  [`inla()`](https://rdrr.io/pkg/INLA/man/inla.html), sometimes after
  altering according to the needs of the inlabru method.

  Warning: Due to how inlabru currently constructs the
  [`inla()`](https://rdrr.io/pkg/INLA/man/inla.html) call, the `mean`,
  `prec`, `mean.intercept`, and `prec.intercept` settings in
  `control.fixed` will have no effect. Until a more elegant alternative
  has been implemented, use explicit `mean.linear` and `prec.linear`
  specifications in each `model="linear"` component instead.

  The following [`inla()`](https://rdrr.io/pkg/INLA/man/inla.html)
  options have inlabru specific defaults:

  `E`

  :   Default `1`.

  `Ntrials`

  :   Default `1L`.

  `control.compute`

  :   Default `list(config = TRUE, control.gcpo = list())`.

  `control.inla`

  :   Default `list(int.strategy = "auto")`.

  `control.fixed`

  :   Default `list(expand.factor.strategy = "inla")`.

## See also

`bru_options()`, `bru_options_default()`, `bru_options_get()`

## Examples

``` r
if (FALSE) { # \dontrun{
if (interactive()) {
  # Combine global and user options:
  options1 <- bru_options(bru_options_get(), bru_verbose = TRUE)
  # Create a proto-options object in two equivalent ways:
  options2 <- as.bru_options(bru_verbose = TRUE)
  options2 <- as.bru_options(list(bru_verbose = TRUE))
  # Combine options objects:
  options3 <- bru_options(options1, options2)
}
} # }
if (FALSE) { # \dontrun{
if (interactive()) {
  bru_options_check(bru_options(bru_max_iter = "text"))
}
} # }
bru_options_get("bru_verbose")
#> [1] 0
if (FALSE) { # \dontrun{
if (interactive()) {
  bru_options_set(
    bru_verbose = TRUE,
    verbose = TRUE
  )
}
} # }
my_fun <- function(val) {
  bru_options_set_local(bru_verbose = val)
  bru_options_get("bru_verbose")
}
# Inside the function, the bru_verbose option is changed.
# Outside the function, the bru_verbose option is unchanged.
print(my_fun(TRUE))
#> [1] TRUE
print(bru_options_get("bru_verbose"))
#> [1] 0
print(my_fun(FALSE))
#> [1] FALSE
print(bru_options_get("bru_verbose"))
#> [1] 0
```
