# Extract standardised names from a bru or inla result object

Extracts the names of fixed effects, random effects, and
hyperparameters, converted with
[`bru_standardise_names()`](https://inlabru-org.github.io/inlabru/reference/bru_standardise_names.md)

## Usage

``` r
bru_names(x)

# S3 method for class 'inla'
bru_names(x)

# S3 method for class 'bru'
bru_names(x)
```

## Arguments

- x:

  an `inla` or
  [`bru()`](https://inlabru-org.github.io/inlabru/reference/bru.md)
  result object

## Value

A character vector with standardised names

## Examples

``` r
if (bru_safe_inla()) {
  fit <- bru(y ~ 1 + x + z(z, model = "iid"),
    data = data.frame(
      y = rnorm(10),
      x = rnorm(10),
      z = rep(seq_len(2), 5)
    )
  )
  bru_names(fit)
}
#>                                         x 
#>                                       "x" 
#>                                 Intercept 
#>                               "Intercept" 
#>                                         z 
#>                                       "z" 
#>   Precision for the Gaussian observations 
#> "Precision_for_the_Gaussian_observations" 
#>                           Precision for z 
#>                         "Precision_for_z" 
```
