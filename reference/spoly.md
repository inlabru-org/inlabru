# Convert a data.frame of boundary points into a SpatialPolgonsDataFrame

A polygon can be described as a sequence of points defining the
polygon's boundary. When given such a sequence (anti clockwise!) this
function creates a `SpatialPolygonsDataFrame` or `sf` holding the
polygon decribed. By default, the first two columns of `data` are
assumed to define the x and y coordinates of the points. This behaviour
can be changed using the `cols` parameter, which points out the names of
the columns holding the coordinates. The coordinate reference system of
the resulting spatial polygon can be set via the `crs` parameter.
Posterior conversion to a different CRS is supported using the `to.crs`
parameter.

## Usage

``` r
spoly(
  data,
  cols = colnames(data)[1:2],
  crs = fm_crs(),
  to.crs = NULL,
  format = c("sp", "sf")
)
```

## Arguments

- data:

  A data.frame of points describing the boundary of the polygon (unique
  points, no holes)

- cols:

  Column names of the x and y coordinates within the data

- crs:

  Coordinate reference system of the points

- to.crs:

  Coordinate reference system for the `SpatialLines`/`sf` ouput.

- format:

  Format of the output object. Either "sp" (default) or "sf"

## Value

[`sp::SpatialPolygonsDataFrame`](https://edzer.github.io/sp/reference/SpatialPolygons.html)
or [sf::sf](https://r-spatial.github.io/sf/reference/sf.html)

## Examples

``` r
# \donttest{
# Create data frame of boundary points (anti clockwise!)
pts <- data.frame(
  x = c(1, 2, 1.7, 1.3),
  y = c(1, 1, 2, 2)
)

# Convert to sf
pol <- spoly(pts, format = "sf")

if (require(ggplot2, quietly = TRUE)) {
  # Plot it!
  ggplot() +
    gg(pol)
}

# }
```
