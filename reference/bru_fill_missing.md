# Fill in missing values in Spatial grids

Computes nearest-available-value imputation for missing values in space

## Usage

``` r
bru_fill_missing(
  data,
  where,
  values,
  layer = NULL,
  selector = NULL,
  batch_size = deprecated()
)
```

## Arguments

- data:

  A SpatialPointsDataFrame, SpatialPixelsDataFrame,
  SpatialGridDataFrame, SpatRaster, Raster, or sf object containing data
  to use for filling

- where:

  A, matrix, data.frame, or SpatialPoints or SpatialPointsDataFrame, or
  sf object, containing the locations of the evaluated values

- values:

  A vector of values to be filled in where `is.na(values)` is `TRUE`

- layer, selector:

  Specifies what data column or columns from which to extract data, see
  [`bru_comp()`](https://inlabru-org.github.io/inlabru/reference/bru_comp.md)
  for details.

- batch_size:

  **\[deprecated\]** due to improved algorithm. Size of
  nearest-neighbour calculation blocks, to limit the memory and
  computational complexity.

## Value

An infilled vector of values

## Examples

``` r
if (FALSE) { # \dontrun{
if (require("sf", quietly = TRUE)) {
  points <-
    sf::st_as_sf(
      data.frame(
        x = 1:3,
        y = 4:6,
        val = c(NA, NA, NA)
      ),
      coords = c("x", "y")
    )
  input_coord <- expand.grid(x = 0:7, y = 0:7)
  input <-
    sf::st_as_sf(
      cbind(input_coord, val = as.vector(input_coord$y)),
      coords = c("x", "y")
    )
  points$val <- bru_fill_missing(input, points, points$val)
  print(points)

  # To fill in missing values in a grid:
  print(input$val[c(3, 30)])
  input$val[c(3, 30)] <- NA # Introduce missing values
  input$val <- bru_fill_missing(input, input, input$val)
  print(input$val[c(3, 30)])
}
} # }
```
