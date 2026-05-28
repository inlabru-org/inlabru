# Grid partitioning

**\[experimental\]** Both the interface and function name may change in
future versions.

Partitions the region based on the given criteria for calculating
residuals in each partition. Parts of this function are taken from
concepts in https://rpubs.com/huanfaChen/grid_from_polygon

## Usage

``` r
cv_partition(
  samplers,
  resolution = NULL,
  nrows = NULL,
  ncols = NULL,
  chess = TRUE,
  ...
)
```

## Arguments

- samplers:

  A sf object containing region for which partitions to be created

- resolution:

  resolution of the grids that are required

- nrows:

  number of rows of grids that are required

- ncols:

  number of columns of grids that are required

- chess:

  chessboard partitioning

- ...:

  Unused.

## Value

a partitioned sf object as required or a list of partitioned sf objects
if chess is TRUE

## Author

Man Ho Suen
