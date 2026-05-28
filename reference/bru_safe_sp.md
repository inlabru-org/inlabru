# Check for potential `sp` version compatibility issues

Loads the sp package with
[`requireNamespace("sp", quietly = TRUE)`](https://rdrr.io/r/base/ns-load.html),
and checks and optionally sets the `sp` evolution status flag if `rgdal`
is unavailable.

## Usage

``` r
bru_safe_sp(quietly = FALSE, force = FALSE, minimum_version = "2.1")
```

## Arguments

- quietly:

  logical; if `TRUE`, prints diagnostic messages. Default `FALSE`

- force:

  logical; If `rgdal` is unavailable and evolution status is less that
  `2L`, return `FALSE` if `force` is `FALSE`. If `force` is `TRUE`,
  return `TRUE` if the package configuration is safe, potentially after
  forcing the evolution status to `2L`. Default `FALSE`

- minimum_version:

  character; the minimum required version. Default 2.1 (should always
  match the requirement in the package DESCRIPTION)

## Value

Returns (invisibly) `FALSE` if a potential issue is detected, and give a
message if `quietly` is `FALSE`. Otherwise returns `TRUE`

## Examples

``` r
if (FALSE) { # \dontrun{
if (bru_safe_sp()) {
  # Run sp dependent calculations
}
} # }
```
