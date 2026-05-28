# Unit test helpers

Local helper functions for package unit tests

## Usage

``` r
local_bru_testthat_assign(x, values, envir = parent.frame())

local_bru_testthat_tolerances(
  tolerances = c(1e-04, 0.01, 0.1),
  envir = parent.frame()
)

local_bru_options_set(..., .reset = FALSE, envir = parent.frame())

local_basic_intercept_testdata()

local_basic_fixed_effect_testdata()

local_bru_safe_inla(multicore = FALSE, quietly = TRUE, envir = parent.frame())

local_bru_testthat_setup(envir = parent.frame())
```

## Arguments

- x:

  character; Name of variable to assign to

- values:

  the object to assign to `x`

- envir:

  environment for exit handlers

- tolerances:

  numeric vector of length 3; `[lowtol, midtol, hitol]`

- .reset:

  For `local_bru_options_set`, logical indicating if the global override
  options list should be emptied before setting the new option(s).

- multicore:

  logical; if `TRUE`, multiple cores are allowed, and the INLA
  `num.threads` option is not checked or altered. Default: `FALSE`,
  multicore not allowed (used for examples and unit tests).

- quietly:

  logical; if `TRUE`, prints diagnostic messages. A message is always
  printed if the INLA `num.threads` option is altered, regardless of the
  `quietly` argument. Default: TRUE.

## Value

`local_bru_options_set()` returns a copy of the global override options
(not including the defaults), invisibly.

## Functions

- `local_bru_testthat_assign()`: Assign local variable. Useful for easy
  cleanup of global workspace with
  [`withr::deferred_run()`](https://withr.r-lib.org/reference/defer.html)
  when running tests interactively.

- `local_bru_testthat_tolerances()`: Assign test tolerances Assign local
  tolerance variables. Useful for easy cleanup of global workspace with
  [`withr::deferred_run()`](https://withr.r-lib.org/reference/defer.html)
  when running tests interactively.

- `local_bru_options_set()`: Wrapper for
  [`bru_options_set_local()`](https://inlabru-org.github.io/inlabru/reference/bru_options.md),
  to locally override the global package options.

- `local_bru_safe_inla()`: Tests should set num.threads = "1:1:1" to
  ensure within-system repeatability by calling `local_bru_safe_inla()`;
  see also
  [`bru_safe_inla()`](https://inlabru-org.github.io/inlabru/reference/bru_safe_inla.md)

- `local_bru_testthat_setup()`: Initialise environment for tests.
  Assigns tolerance variables. To be called either at the top of a
  testfile, or inside tests. Does *not* call `local_bru_safe_inla()`,
  since that may invoke a skip and should be called inside each test
  that relies on INLA.

## See also

[`bru_options_set_local()`](https://inlabru-org.github.io/inlabru/reference/bru_options.md),
[`bru_options_default()`](https://inlabru-org.github.io/inlabru/reference/bru_options.md),
[`bru_options_get()`](https://inlabru-org.github.io/inlabru/reference/bru_options.md)
