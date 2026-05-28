# Geom for fm_mesh_2d objects

This function extracts the graph of an
[fmesher::fm_mesh_2d](https://inlabru-org.github.io/fmesher/reference/fm_mesh_2d.html)
object and uses `geom_line` to visualize the graph's edges.
Alternatively, if the `color` argument is provided, interpolates the
colors across for a set of SpatialPixels covering the mesh area and
calls
[`gg.SpatialPixelsDataFrame()`](https://inlabru-org.github.io/inlabru/reference/gg.Spatial.md)
to plot the interpolation. Requires the `ggplot2` package.

Also see the
[`fmesher::geom_fm()`](https://inlabru-org.github.io/fmesher/reference/geom_fm.html)
method.

## Usage

``` r
# S3 method for class 'fm_mesh_2d'
gg(
  data,
  color = NULL,
  alpha = NULL,
  edge.color = "grey",
  edge.linewidth = 0.25,
  interior = TRUE,
  int.color = "blue",
  int.linewidth = 0.5,
  exterior = TRUE,
  ext.color = "black",
  ext.linewidth = 1,
  crs = NULL,
  mask = NULL,
  nx = 500,
  ny = 500,
  ...
)
```

## Arguments

- data:

  An `fm_mesh_2d` object.

- color:

  A vector of scalar values to fill the mesh with colors. The length of
  the vector mus correspond to the number of mesh vertices. The
  alternative name `colour` is also recognised.

- alpha:

  A vector of scalar values setting the alpha value of the colors
  provided.

- edge.color:

  Color of the regular mesh edges.

- edge.linewidth:

  Line width for the regular mesh edges. Default 0.25

- interior:

  If TRUE, plot the interior boundaries of the mesh.

- int.color:

  Color used to plot the interior constraint edges.

- int.linewidth:

  Line width for the interior constraint edges. Default 0.5

- exterior:

  If TRUE, plot the exterior boundaries of the mesh.

- ext.color:

  Color used to plot the exterior boundary edges.

- ext.linewidth:

  Line width for the exterior boundary edges. Default 1

- crs:

  A CRS object supported by
  [`fmesher::fm_transform()`](https://inlabru-org.github.io/fmesher/reference/fm_transform.html)
  defining the coordinate system to project the mesh to before plotting.

- mask:

  A `SpatialPolygon` or `sf` polygon defining the region that is
  plotted.

- nx:

  Number of pixels in x direction (when plotting using the color
  parameter).

- ny:

  Number of pixels in y direction (when plotting using the color
  parameter).

- ...:

  ignored arguments (S3 generic compatibility).

## Value

`geom_line` return values or, if the color argument is used, the values
of
[`gg.SpatialPixelsDataFrame()`](https://inlabru-org.github.io/inlabru/reference/gg.Spatial.md).

## See also

Other geomes:
[`gg()`](https://inlabru-org.github.io/inlabru/reference/gg.md),
[`gg.RasterLayer()`](https://inlabru-org.github.io/inlabru/reference/gg.RasterLayer.md),
[`gg.SpatRaster()`](https://inlabru-org.github.io/inlabru/reference/gg.SpatRaster.md),
[`gg.Spatial`](https://inlabru-org.github.io/inlabru/reference/gg.Spatial.md),
[`gg.data.frame()`](https://inlabru-org.github.io/inlabru/reference/gg.bru_prediction.md),
[`gg.fm_mesh_1d()`](https://inlabru-org.github.io/inlabru/reference/gg.fm_mesh_1d.md),
[`gg.matrix()`](https://inlabru-org.github.io/inlabru/reference/gg.matrix.md),
[`gg.sf()`](https://inlabru-org.github.io/inlabru/reference/gg.sf.md)

## Examples

``` r
# \donttest{
if (require(fmesher, quietly = TRUE) &&
  require(ggplot2, quietly = TRUE)) {
  # Load Gorilla data
  gorillas <- inlabru::gorillas_sf

  # Plot mesh using default edge colors

  ggplot() +
    gg(gorillas$mesh)

  # Don't show interior and exterior boundaries

  ggplot() +
    gg(gorillas$mesh, interior = FALSE, exterior = FALSE)

  # Change the edge colors

  ggplot() +
    gg(gorillas$mesh,
      edge.color = "green",
      int.color = "black",
      ext.color = "blue"
    )

  # Use the x-coordinate of the vertices to colorize the triangles and
  # mask the plotted area by the survey boundary, i.e. only plot the inside

  xcoord <- gorillas$mesh$loc[, 1]
  ggplot() +
    gg(gorillas$mesh, color = (xcoord - 580), mask = gorillas$boundary) +
    gg(gorillas$boundary, alpha = 0)
}

# }
```
