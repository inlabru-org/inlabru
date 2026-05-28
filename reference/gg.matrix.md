# Geom for matrix

Creates a tile geom for plotting a matrix

## Usage

``` r
# S3 method for class 'matrix'
gg(data, mapping = NULL, ...)
```

## Arguments

- data:

  A `matrix` object.

- mapping:

  a set of aesthetic mappings created by `aes`. These are passed on to
  `geom_tile`.

- ...:

  Arguments passed on to `geom_tile`.

## Value

A `geom_tile` with reversed y scale.

## Details

Requires the `ggplot2` package.

## See also

Other geomes:
[`gg()`](https://inlabru-org.github.io/inlabru/reference/gg.md),
[`gg.RasterLayer()`](https://inlabru-org.github.io/inlabru/reference/gg.RasterLayer.md),
[`gg.SpatRaster()`](https://inlabru-org.github.io/inlabru/reference/gg.SpatRaster.md),
[`gg.Spatial`](https://inlabru-org.github.io/inlabru/reference/gg.Spatial.md),
[`gg.data.frame()`](https://inlabru-org.github.io/inlabru/reference/gg.bru_prediction.md),
[`gg.fm_mesh_1d()`](https://inlabru-org.github.io/inlabru/reference/gg.fm_mesh_1d.md),
[`gg.fm_mesh_2d()`](https://inlabru-org.github.io/inlabru/reference/gg.fm_mesh_2d.md),
[`gg.sf()`](https://inlabru-org.github.io/inlabru/reference/gg.sf.md)

## Examples

``` r
if (require("ggplot2", quietly = TRUE)) {
  A <- matrix(runif(100), nrow = 10)
  ggplot() +
    gg(A)
}
```
