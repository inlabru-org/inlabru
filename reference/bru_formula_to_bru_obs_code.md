# Convert components to R code

Convert a [formula](https://rdrr.io/r/stats/formula.html) describing
latent model components to R code strings that can be evaluated to
create the corresponding
[bru_comp](https://inlabru-org.github.io/inlabru/reference/bru_comp.md)
objects.

## Usage

``` r
bru_formula_to_bru_obs_code(components, add = "")
```

## Arguments

- components:

  A [formula](https://rdrr.io/r/stats/formula.html) describing latent
  model components.

## Value

a character vector of R code strings, one for each component in the
formula.

## Author

Fabian E. Bachl <bachlfab@gmail.com>
