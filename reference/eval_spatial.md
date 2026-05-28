# Evaluate spatial covariates

Evaluate spatial covariates

## Usage

``` r
eval_spatial(data, where, layer = NULL, selector = NULL)

# S3 method for class 'SpatialPolygonsDataFrame'
eval_spatial(data, where, layer = NULL, selector = NULL)

# S3 method for class 'SpatialPixelsDataFrame'
eval_spatial(data, where, layer = NULL, selector = NULL)

# S3 method for class 'SpatialGridDataFrame'
eval_spatial(data, where, layer = NULL, selector = NULL)

# S3 method for class 'sf'
eval_spatial(data, where, layer = NULL, selector = NULL)

# S3 method for class 'SpatRaster'
eval_spatial(data, where, layer = NULL, selector = NULL)

# S3 method for class 'stars'
eval_spatial(data, where, layer = NULL, selector = NULL)
```

## Arguments

- data:

  Spatial data

- where:

  Where to evaluate the data

- layer:

  Which `data` layer to extract (as integer or character). May be a
  vector, specifying a separate layer for each `where` item.

- selector:

  The name of a variable in `where` specifying the `layer` information.

## Methods (by class)

- `eval_spatial(SpatialPolygonsDataFrame)`: Compatibility wrapper for
  `eval_spatial.sf`

- `eval_spatial(sf)`: Supports point-in-polygon information lookup.
  Other combinations are untested.
