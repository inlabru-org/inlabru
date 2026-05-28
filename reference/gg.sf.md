# Geom helper for sf objects

This function uses
[`geom_sf()`](https://ggplot2.tidyverse.org/reference/ggsf.html), unless
overridden by the geom argument. Requires the `ggplot2` package.

## Usage

``` r
# S3 method for class 'sf'
gg(data, mapping = NULL, ..., geom = "sf")
```

## Arguments

- data:

  An `sf` object.

- mapping:

  Default mapping is `ggplot2::aes(geometry = ...)`, where the geometry
  name is obtained from `attr(data, "sf_column")`. This is merged with
  the user supplied mapping.

- ...:

  Arguments passed on to `geom_sf` or `geom_tile`.

- geom:

  Either "sf" (default) or "tile". For "tile", uses
  `geom_tile(..., stat = "sf_coordinates")`, intended for converting
  point data to grid tiles with the `fill` aesthetic, which is by
  default set to the first data column.

## Value

A ggplot return value

## See also

Other geomes:
[`gg()`](https://inlabru-org.github.io/inlabru/reference/gg.md),
[`gg.RasterLayer()`](https://inlabru-org.github.io/inlabru/reference/gg.RasterLayer.md),
[`gg.SpatRaster()`](https://inlabru-org.github.io/inlabru/reference/gg.SpatRaster.md),
[`gg.Spatial`](https://inlabru-org.github.io/inlabru/reference/gg.Spatial.md),
[`gg.data.frame()`](https://inlabru-org.github.io/inlabru/reference/gg.bru_prediction.md),
[`gg.fm_mesh_1d()`](https://inlabru-org.github.io/inlabru/reference/gg.fm_mesh_1d.md),
[`gg.fm_mesh_2d()`](https://inlabru-org.github.io/inlabru/reference/gg.fm_mesh_2d.md),
[`gg.matrix()`](https://inlabru-org.github.io/inlabru/reference/gg.matrix.md)

## Examples

``` r
# \donttest{
if (require("ggplot2", quietly = TRUE) &&
  bru_safe_terra(quietly = TRUE) &&
  require("tidyterra", quietly = TRUE)) {
  # Load Gorilla data

  gorillas <- inlabru::gorillas_sf
  gorillas$gcov <- gorillas_sf_gcov()

  # Plot Gorilla elevation covariate provided as terra::rast.

  ggplot() +
    gg(gorillas$gcov$elevation)

  # Add Gorilla survey boundary and nest sightings

  ggplot() +
    gg(gorillas$gcov$elevation) +
    gg(gorillas$boundary, alpha = 0) +
    gg(gorillas$nests)

  # Load pantropical dolphin data

  mexdolphin <- inlabru::mexdolphin_sf

  # Plot the pantropical survey boundary, ship transects and dolphin
  # sightings

  ggplot() +
    gg(mexdolphin$ppoly, alpha = 0.5) + # survey boundary
    gg(mexdolphin$samplers) + # ship transects
    gg(mexdolphin$points) # dolphin sightings

  # Change color

  ggplot() +
    gg(mexdolphin$ppoly, color = "green", alpha = 0.5) + # survey boundary
    gg(mexdolphin$samplers, color = "red") + # ship transects
    gg(mexdolphin$points, color = "blue") # dolphin sightings


  # Visualize data annotations: line width by segment number

  names(mexdolphin$samplers) # 'seg' holds the segment number
  ggplot() +
    gg(mexdolphin$samplers, aes(color = seg))

  # Visualize data annotations: point size by dolphin group size

  names(mexdolphin$points) # 'size' holds the group size
  ggplot() +
    gg(mexdolphin$points, aes(size = size))
}

# }
```
