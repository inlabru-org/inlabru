# Convert SpatialPoints and boundary polygon to spatstat ppp object

Spatstat point pattern objects consist of points and an observation
windows. This function uses a `SpatialPoints` object and a
`SpatialPolygon` object to generate the points and the window. Lastly,
the `ppp()` function is called to create the `ppp` object.

## Usage

``` r
spatial.to.ppp(points, samplers)
```

## Arguments

- points:

  A `SpatialPoints[DataFrame]` object describing the point pattern.

- samplers:

  A `SpatialPolygons[DataFrame]` object describing the observation
  window.

## Value

A spatstat `spatstat` `ppp` object

## Examples

``` r
# \donttest{
if (require("spatstat.geom") &&
  bru_safe_sp() &&
  require("sp") &&
  bru_safe_terra(quietly = TRUE) &&
  require("sf", quietly = TRUE)) {
  # Load Gorilla data

  gorillas <- gorillas_sp()

  # Use nest locations and survey boundary to create a spatstat ppp object

  gp <- spatial.to.ppp(gorillas$nests, gorillas$boundary)
  class(gp)

  # Plot it

  plot(gp)
}
#> Loading required package: spatstat.geom
#> Loading required package: spatstat.data
#> Loading required package: spatstat.univar
#> spatstat.univar 3.1-7
#> spatstat.geom 3.7-3
#> 
#> Attaching package: ‘spatstat.geom’
#> The following objects are masked from ‘package:raster’:
#> 
#>     area, rotate, shift
#> The following object is masked from ‘package:patchwork’:
#> 
#>     area
#> Warning: data contain duplicated points

# }
```
