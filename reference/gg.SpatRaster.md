# Geom wrapper for SpatRaster objects

Convenience wrapper function for
[`tidyterra::geom_spatraster()`](https://dieghernan.github.io/tidyterra/reference/geom_spatraster.html).
Requires the `ggplot2` and `tidyterra` packages.

## Usage

``` r
# S3 method for class 'SpatRaster'
gg(data, ...)
```

## Arguments

- data:

  A SpatRaster object.

- ...:

  Arguments passed on to `geom_spatraster`.

## Value

The output from \`geom_spatraster.

## See also

Other geomes:
[`gg()`](https://inlabru-org.github.io/inlabru/reference/gg.md),
[`gg.RasterLayer()`](https://inlabru-org.github.io/inlabru/reference/gg.RasterLayer.md),
[`gg.Spatial`](https://inlabru-org.github.io/inlabru/reference/gg.Spatial.md),
[`gg.data.frame()`](https://inlabru-org.github.io/inlabru/reference/gg.bru_prediction.md),
[`gg.fm_mesh_1d()`](https://inlabru-org.github.io/inlabru/reference/gg.fm_mesh_1d.md),
[`gg.fm_mesh_2d()`](https://inlabru-org.github.io/inlabru/reference/gg.fm_mesh_2d.md),
[`gg.matrix()`](https://inlabru-org.github.io/inlabru/reference/gg.matrix.md),
[`gg.sf()`](https://inlabru-org.github.io/inlabru/reference/gg.sf.md)

## Examples

``` r
# \donttest{
if (require("ggplot2", quietly = TRUE) &&
  bru_safe_terra(quietly = TRUE) &&
  require("tidyterra", quietly = TRUE)) {
  # Load Gorilla covariates

  gcov <- gorillas_sf_gcov()

  # Plot the pixel centers
  ggplot() +
    gg(gcov$elevation)
}
#> 
#> Attaching package: ‘tidyterra’
#> The following object is masked from ‘package:stats’:
#> 
#>     filter

# }
```
