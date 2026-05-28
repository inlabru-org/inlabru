# Load INLA safely for examples and tests

Loads the INLA package with
[`requireNamespace("INLA", quietly = TRUE)`](https://rdrr.io/r/base/ns-load.html),
and optionally checks and sets the multicore `num.threads` INLA option.

## Usage

``` r
bru_safe_inla(multicore = NULL, quietly = FALSE, minimum_version = "23.1.31")
```

## Arguments

- multicore:

  logical; if `TRUE`, multiple cores are allowed, and the INLA
  `num.threads` option is not checked or altered. If `FALSE`, forces
  `num.threads="1:1:1"`. Default: NULL, checks if running in testthat or
  non-interactively, in which case sets `multicore=FALSE`, otherwise
  `TRUE`.

- quietly:

  logical; if `FALSE` and `multicore` is `FALSE`, prints a message if
  the `num.threads` option isn't already "1.1:1" to alert the user to
  the change. Default: FALSE.

- minimum_version:

  character; the minimum required INLA version. Default 23.1.31 (should
  always match the requirement in the package DESCRIPTION)

## Value

logical; `TRUE` if INLA was loaded safely, otherwise FALSE

## Examples

``` r
if (FALSE) { # \dontrun{
if (bru_safe_inla()) {
  # Run inla dependent calculations
}
} # }
```
