# Geom for RasterLayer objects

This function takes a RasterLayer object, converts it into a
`SpatialPixelsDataFrame` and uses `geom_tile` to plot the data.

## Usage

``` r
# S3 method for class 'RasterLayer'
gg(
  data,
  mapping = ggplot2::aes(x = .data[["x"]], y = .data[["y"]], fill = .data[["layer"]]),
  ...
)
```

## Arguments

- data:

  A RasterLayer object.

- mapping:

  aesthetic mappings created by `aes`. These are passed on to
  `geom_tile`.

- ...:

  Arguments passed on to `geom_tile`.

## Value

An object returned by `geom_tile`

## Details

This function requires the `raster` and `ggplot2` packages.

## See also

Other geomes:
[`gg()`](https://inlabru-org.github.io/inlabru/reference/gg.md),
[`gg.SpatRaster()`](https://inlabru-org.github.io/inlabru/reference/gg.SpatRaster.md),
[`gg.Spatial`](https://inlabru-org.github.io/inlabru/reference/gg.Spatial.md),
[`gg.data.frame()`](https://inlabru-org.github.io/inlabru/reference/gg.bru_prediction.md),
[`gg.fm_mesh_1d()`](https://inlabru-org.github.io/inlabru/reference/gg.fm_mesh_1d.md),
[`gg.fm_mesh_2d()`](https://inlabru-org.github.io/inlabru/reference/gg.fm_mesh_2d.md),
[`gg.matrix()`](https://inlabru-org.github.io/inlabru/reference/gg.matrix.md),
[`gg.sf()`](https://inlabru-org.github.io/inlabru/reference/gg.sf.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Some features require the raster and spatstat.data packages.
if (require("spatstat.data", quietly = TRUE) &&
  require("raster", quietly = TRUE) &&
  require("ggplot2", quietly = TRUE)) {
  # Load Gorilla data
  data("gorillas", package = "spatstat.data", envir = environment())

  # Convert elevation covariate to RasterLayer

  elev <- as(gorillas.extra$elevation, "RasterLayer")

  # Plot the elevation

  ggplot() +
    gg(elev)
}
} # }
```
