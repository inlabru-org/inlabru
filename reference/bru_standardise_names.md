# Standardise inla hyperparameter names

The inla hyperparameter output uses parameter names that can include
whitespace and special characters. This function replaces those
characters with underscores.

## Usage

``` r
bru_standardise_names(x)
```

## Arguments

- x:

  character vector; names to be standardised

## Value

A character vector with standardised names

## See also

[`bru_names()`](https://inlabru-org.github.io/inlabru/reference/bru_names.md)

## Examples

``` r
bru_standardise_names("Precision for the Gaussian observations")
#>   Precision for the Gaussian observations 
#> "Precision_for_the_Gaussian_observations" 
```
