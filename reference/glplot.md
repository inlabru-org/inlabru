# Render objects using RGL

`glplot()` is a generic function for renders various kinds of spatial
objects, i.e. `Spatial*` data and `fm_mesh_2d` objects. The function
invokes particular methods which depend on the class of the first
argument.

## Usage

``` r
glplot(object, ...)

# S3 method for class 'SpatialPoints'
glplot(object, add = TRUE, color = "red", ...)

# S3 method for class 'SpatialLines'
glplot(object, add = TRUE, ...)

# S3 method for class 'fm_mesh_2d'
glplot(object, add = TRUE, col = NULL, ...)

globe(
  R = 1,
  R.grid = 1.05,
  specular = "black",
  axes = FALSE,
  box = FALSE,
  xlab = "",
  ylab = "",
  zlab = ""
)
```

## Arguments

- object:

  an object used to select a method.

- ...:

  Parameters passed on to `plot_rgl.fm_mesh_2d()`

- add:

  If TRUE, add the points to an existing plot. If FALSE, create new
  plot.

- color:

  vector of R color characters. See material3d() for details.

- col:

  Color specification. A single named color, a vector of scalar values,
  or a matrix of RGB values.

- R:

  Radius of the globe

- R.grid:

  Radius of the annotation sphere.

- specular:

  Light color of specular effect.

- axes:

  If TRUE, plot x, y and z axes.

- box:

  If TRUE, plot a box around the globe.

- xlab, ylab, zlab:

  Axes labels

## Methods (by class)

- `glplot(SpatialPoints)`: This function will calculate the cartesian
  coordinates of the points provided and use points3d() in order to
  render them.

- `glplot(SpatialLines)`: This function will calculate a cartesian
  representation of the lines provided and use
  [`lines3d()`](https://dmurdoch.github.io/rgl/dev/reference/primitives.html)
  in order to render them.

- `glplot(fm_mesh_2d)`: This function transforms the mesh to 3D
  cartesian coordinates and uses
  [`fmesher::plot_rgl()`](https://inlabru-org.github.io/fmesher/reference/plot_rgl.html)
  to plot the result.

## Functions

- `globe()`: Visualize a globe using RGL

  Creates a textured sphere and lon/lat coordinate annotations. This
  function requires the `rgl` and `sphereplot` packages.

## Examples

``` r
# \donttest{
if (interactive() &&
  require("rgl", quietly = TRUE) &&
  require("sphereplot", quietly = TRUE) &&
  bru_safe_sp() &&
  require("sp")) {
  # Show the globe:
  globe()

  # Load pantropoical dolphin data
  mexdolphin <- inlabru::mexdolphin_sp()

  # Add mesh, ship transects and dolphin sightings stored
  # as fm_mesh_2d, SpatialLines and SpatialPoints objects, respectively

  glplot(mexdolphin$mesh, alpha = 0.2)
  glplot(mexdolphin$samplers, lwd = 5)
  glplot(mexdolphin$points, size = 10)
}
# }
```
