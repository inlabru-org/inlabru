# Add component input/latent mappers

Add missing mappers between input data and latent variables, based on
likelihood data

Equip component(s) with mappers for subcomponents that do not have
predefined mappers. When needed, the data in `lhoods` is used to
determine the appropriate mapper(s).

## Usage

``` r
add_mappers(...)

# S3 method for class 'bru_comp'
add_mappers(component, lhoods, inputs = NULL, lh_data = NULL, ...)

# S3 method for class 'bru_comp_list'
add_mappers(components, lhoods, inputs = NULL, lh_data = NULL, ...)
```

## Arguments

- ...:

  Parameters passed on to other methods

- component:

  A
  [component](https://inlabru-org.github.io/inlabru/reference/bru_comp.md)
  object

- lhoods:

  A
  [bru_obs_list](https://inlabru-org.github.io/inlabru/reference/bru_obs.md)
  object

- inputs:

  A precomputed list of inputs, as returned by
  [`bru_input.bru_comp_list()`](https://inlabru-org.github.io/inlabru/reference/bru_input.md)

- lh_data:

  A list of data object, one for each likelihood in `lhoods` that is
  used to determine the mapper(s). If `NULL`, the data is or inputs are
  used instead.

- components:

  A `bru_comp_list` object

## Value

A `component` object with completed mapper information

## Examples

``` r
if (FALSE) { # \dontrun{
if (interactive() && bru_safe_inla()) {}
} # }
```
