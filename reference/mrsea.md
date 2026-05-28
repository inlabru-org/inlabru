# Marine renewables strategic environmental assessment

Data imported from package MRSea, see
<https://www.creem.st-andrews.ac.uk/software/>

## Usage

``` r
mrsea
```

## Format

A list of objects:

- `points`:

  A `sf` object containing the locations of XXXXX.

- `samplers`:

  A `sf` object containing the transect lines that were surveyed.

- `mesh`:

  An `fm_mesh_2d` object containing a Delaunay triangulation mesh (a
  type of discretization of continuous space) covering the survey
  region.

- `boundary`:

  An `sf` object defining the boundary polygon of the survey region.

- `covar`:

  An `sf` containing sea depth estimates.

## Source

Library `MRSea`.

## References

NONE YET

## Examples

``` r
if (require(ggplot2, quietly = TRUE)) {
  ggplot() +
    fmesher::geom_fm(data = mrsea$mesh) +
    gg(mrsea$samplers) +
    gg(mrsea$points) +
    gg(mrsea$boundary)
}
```
