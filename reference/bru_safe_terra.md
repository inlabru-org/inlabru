# Check for potential `terra` installation issues

Loads the terra package with
[`requireNamespace("terra", quietly = TRUE)`](https://rdrr.io/r/base/ns-load.html),
and checks whether
[`eval_spatial()`](https://inlabru-org.github.io/inlabru/reference/eval_spatial.md)
successfully extracts non-NA values. NA values may appear when there is
a GDAL/PROJ installation issue that causes it to not find the `proj.db`
database file.

For use in package tests and examples that depend on `terra`.

## Usage

``` r
bru_safe_terra(quietly = FALSE, minimum_version = "1.7-66")
```

## Arguments

- quietly:

  logical; if `TRUE`, prints diagnostic messages. Default `FALSE`

- minimum_version:

  character; the minimum required version. Default 1.7-66 (should always
  match the requirement in the package DESCRIPTION)

## Value

Returns (invisibly) `FALSE` if a potential issue is detected, and give a
message if `quietly` is `FALSE`. Otherwise returns `TRUE`

## Examples

``` r
if (FALSE) { # \dontrun{
if (bru_safe_terra()) {
  # Run terra dependent calculations
}
} # }
```
