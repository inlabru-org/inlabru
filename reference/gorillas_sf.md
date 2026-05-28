# Gorilla nesting sites in sf format

This is the `gorillas` dataset from the package `spatstat.data`,
reformatted as point process data for use with `inlabru`.

## Usage

``` r
gorillas_sf
data(gorillas_sf, package = "inlabru")

gorillas_sf_gcov()

gorillas_sp()
```

## Format

The data are a list that contains these elements:

- `nests`::

  An `sf` object containing the locations of the gorilla nests.

- `boundary`::

  An `sf` object defining the boundary of the region that was searched
  for the nests.

- `mesh`::

  An `fm_mesh_2d` object containing a mesh that can be used with
  function `lgcp` to fit a LGCP to the nest data.

- `gcov_file`::

  The in-package filename of a
  [`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  object, with one layer for each of these spatial covariates:

  `aspect`

  :   Compass direction of the terrain slope. Categorical, with levels
      N, NE, E, SE, S, SW, W and NW, which are coded as integers 1 to 8.

  `elevation`

  :   Digital elevation of terrain, in metres.

  `heat`

  :   Heat Load Index at each point on the surface (Beer's aspect),
      discretised. Categorical with values Warmest (Beer's aspect
      between 0 and 0.999), Moderate (Beer's aspect between 1 and
      1.999), Coolest (Beer's aspect equals 2). These are coded as
      integers 1, 2 and 3, in that order.

  `slopangle`

  :   Terrain slope, in degrees.

  `slopetype`

  :   Type of slope. Categorical, with values Valley, Toe (toe slope),
      Flat, Midslope, Upper and Ridge. These are coded as integers 1 to
      6.

  `vegetation`

  :   Vegetation type: a categorical variable with 6 levels coded as
      integers 1 to 6 (in order of increasing expected habitat
      suitability)

  `waterdist`

  :   Euclidean distance from nearest water body, in metres.

  Loading of the covariates can be done with `gorillas_sf_gcov()` or

      gorillas_sf$gcov <- terra::rast(
        system.file(gorillas_sf$gcov_file, package = "inlabru")
      )

- `plotsample`:

  Plot sample of gorilla nests, sampling 9x9 over the region, with 60\\

  `counts`

  :   An `sf` object with elements `count`, `exposure`, and `geometry`,
      holding the point geometry for the centre of each plot, the count
      in each plot and the area of each plot.

  `plots`

  :   An `sf` object with `MULTIPOLYGON` objects defining the individual
      plot boundaries and an all-ones `weight` column.

  `nests`

  :   An `sf` giving the locations of each detected nests, `group`
      ("minor" or "major"), `season` ("dry" or "rainy"), and `date` (in
      `Date` format).

## Source

Library `spatstat.data`.

## Functions

- `gorillas_sf_gcov()`: Access the `gorillas_sf` covariates data as a
  [`terra::rast()`](https://rspatial.github.io/terra/reference/rast.html)
  object.

- `gorillas_sp()`: Access the `gorillas_sf` data in `sp` format. The
  covariate data is added as `gcov`, a list of
  [`sp::SpatialPixelsDataFrame`](https://edzer.github.io/sp/reference/SpatialGridDataFrame.html)
  objects. Requires the `sp`, `sf`, and `terra` packages to be
  installed.

## References

Funwi-Gabga, N. (2008) A pastoralist survey and fire impact assessment
in the Kagwene Gorilla Sanctuary, Cameroon. M.Sc. thesis, Geology and
Environmental Science, University of Buea, Cameroon.

Funwi-Gabga, N. and Mateu, J. (2012) Understanding the nesting spatial
behaviour of gorillas in the Kagwene Sanctuary, Cameroon. Stochastic
Environmental Research and Risk Assessment 26 (6), 793-811.

## Examples

``` r
if (interactive() &&
  require(ggplot2, quietly = TRUE) &&
  bru_safe_terra(quietly = TRUE) &&
  requireNamespace("tidyterra", quietly = TRUE)) {
  # plot all the nests, mesh and boundary
  ggplot() +
    gg(gorillas_sf$mesh) +
    geom_sf(
      data = gorillas_sf$boundary,
      alpha = 0.1, fill = "blue"
    ) +
    geom_sf(data = gorillas_sf$nests)

  # Plot the elevation covariate
  gorillas_sf$gcov <- gorillas_sf_gcov()
  ggplot() +
    tidyterra::geom_spatraster(data = gorillas_sf$gcov$elevation)

  # Plot the plot sample
  ggplot() +
    geom_sf(data = gorillas_sf$plotsample$plots) +
    geom_sf(data = gorillas_sf$plotsample$nests)
}
if (interactive() &&
  bru_safe_terra(quietly = TRUE)) {
  gorillas_sf$gcov <- gorillas_sf_gcov()
}
```
