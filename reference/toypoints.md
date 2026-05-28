# Simulated 2D point process data

This data set serves as an example for basic inlabru.

## Usage

``` r
data(toypoints)
```

## Format

The data are a list that contains these elements:

- `points`:

  An `sf` object of point locations and and `z` measurements

- `mesh`:

  An `fm_mesh_2d` object

- `boundary`:

  An `sf` polygon denting the region of interest

- `pred_locs`:

  A `sf` object with prediction point locations

## Examples

``` r
if (require("ggplot2")) {
  ggplot() +
    fmesher::geom_fm(data = toypoints$mesh, alpha = 0) +
    geom_sf(data = toypoints$boundary, fill = "blue", alpha = 0.1) +
    geom_sf(data = toypoints$points, aes(color = z)) +
    scale_color_viridis_c()
}
```
