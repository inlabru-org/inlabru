# Blue and red shrimp in the Western Mediterranean Sea

Blue and red shrimp in the Western Mediterranean Sea.

## Usage

``` r
data(shrimp)
```

## Format

A list of objects:

- `hauls`::

  An `sf` object containing haul locations

- `mesh`::

  An `fm_mesh_2d` object containing a Delaunay triangulation mesh (a
  type of discretization of continuous space) covering the haul
  locations.

## Source

Pennino, Maria Grazia. Personal communication.

## References

Pennino, M. G., Paradinas, I., Munoz, F., Illian, J.,Quilez-Lopez, A.,
Bellido, J.M., Conesa, D. Accounting for preferential sampling in
species distribution models. Ecology and Evolution, 9(1), p653-663, 2019
[doi:10.1002/ece3.4789](https://doi.org/10.1002/ece3.4789)

## Examples

``` r
# \donttest{
if (require(ggplot2, quietly = TRUE)) {
  data(shrimp, package = "inlabru", envir = environment())
  ggplot() +
    fmesher::geom_fm(data = shrimp$mesh) +
    gg(shrimp$hauls, aes(col = catch)) +
    coord_sf(datum = fmesher::fm_crs(shrimp$hauls))
}

# }
```
