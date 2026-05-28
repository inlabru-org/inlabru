# ggplot2 geomes for inlabru related objects

gg is a generic function for generating geomes from various kinds of
spatial objects, e.g. Spatial\* data, meshes, Raster objects and
inla/inlabru predictions. The function invokes particular methods which
depend on the [class](https://rdrr.io/r/base/class.html) of the first
argument.

## Usage

``` r
gg(data, ...)
```

## Arguments

- data:

  an object for which to generate a geom.

- ...:

  Arguments passed on to the geom method.

## Value

The form of the value returned by gg depends on the class of its
argument. See the documentation of the particular methods for details of
what is produced by that method.

## See also

Other geomes:
[`gg.RasterLayer()`](https://inlabru-org.github.io/inlabru/reference/gg.RasterLayer.md),
[`gg.SpatRaster()`](https://inlabru-org.github.io/inlabru/reference/gg.SpatRaster.md),
[`gg.Spatial`](https://inlabru-org.github.io/inlabru/reference/gg.Spatial.md),
[`gg.data.frame()`](https://inlabru-org.github.io/inlabru/reference/gg.bru_prediction.md),
[`gg.fm_mesh_1d()`](https://inlabru-org.github.io/inlabru/reference/gg.fm_mesh_1d.md),
[`gg.fm_mesh_2d()`](https://inlabru-org.github.io/inlabru/reference/gg.fm_mesh_2d.md),
[`gg.matrix()`](https://inlabru-org.github.io/inlabru/reference/gg.matrix.md),
[`gg.sf()`](https://inlabru-org.github.io/inlabru/reference/gg.sf.md)

## Examples

``` r
if (require("ggplot2", quietly = TRUE)) {
  # Load Gorilla data

  gorillas <- inlabru::gorillas_sf

  # Invoke ggplot and add geomes for the Gorilla nests and the survey
  # boundary

  ggplot() +
    gg(gorillas$boundary) +
    gg(gorillas$nests)
}
```
