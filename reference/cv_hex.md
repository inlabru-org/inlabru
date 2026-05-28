# Hexagon tiling of region

**\[experimental\]** Both the interface and function name may change in
future versions.

Partitions the region based on the hexagon tiling

Input:

## Usage

``` r
cv_hex(samplers, cellsize = 0.5, n_group = 3, ...)
```

## Arguments

- samplers:

  A sf object containing region for which partitions to be created

- cellsize:

  hexagon cellsize, see
  [`sf::st_make_grid()`](https://r-spatial.github.io/sf/reference/st_make_grid.html)
  description

- n_group:

  number of cv folds.

- ...:

  Passed on to
  [`fm_nonconvex_hull()`](https://inlabru-org.github.io/fmesher/reference/fm_nonconvex_hull.html),
  e.g. `resolution`

## Value

a hexagonal partition covering the `samplers` object, with a `group`
integer variable indicating the fold assignment, 1, 2, ..., `n_group`.

## Author

Man Ho Suen

## Examples

``` r
if (interactive()) {
  bnd <- gorillas_sf$boundary
  hex_cv <- cv_hex(bnd, cellsize = 0.5, n_group = 3, resolution = 100)
  plot(hex_cv)

  chess <- cv_partition(bnd, resolution = 0.5, chess = TRUE)
  plot(chess$white)
}
```
