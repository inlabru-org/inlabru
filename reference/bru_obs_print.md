# Summary and print methods for observation models

Summary and print methods for observation models

## Usage

``` r
# S3 method for class 'bru_obs'
summary(object, verbose = TRUE, ...)

# S3 method for class 'bru_obs_list'
summary(object, verbose = TRUE, ...)

# S3 method for class 'summary_bru_obs'
print(x, ...)

# S3 method for class 'summary_bru_obs_list'
print(x, ...)

# S3 method for class 'bru_obs'
print(x, ...)

# S3 method for class 'bru_obs_list'
print(x, ...)
```

## Arguments

- object:

  Object to operate on

- verbose:

  logical; If `TRUE`, include more details of the component definitions.
  If `FALSE`, only show basic component definition information. Default:
  `TRUE`

- ...:

  Arguments passed on to other `summary` methods

- x:

  Object to be printed

## See also

[`bru_obs()`](https://inlabru-org.github.io/inlabru/reference/bru_obs.md)

## Examples

``` r
obs <- bru_obs(y ~ ., data = data.frame(y = rnorm(10)))
summary(obs)
#>   Model tag: <No tag>
#>     Family: 'gaussian'
#>     Data class: 'data.frame'
#>     Response class: 'numeric'
#>     Predictor: y ~ .
#>     Additive/Linear/Rowwise: TRUE/TRUE/TRUE
#>     Used components: effect[<not yet initialised>], latent[<not yet initialised>] 
print(obs)
#>   Model tag: <No tag>
#>     Family: 'gaussian'
#>     Data class: 'data.frame'
#>     Response class: 'numeric'
#>     Predictor: y ~ .
#>     Additive/Linear/Rowwise: TRUE/TRUE/TRUE
#>     Used components: effect[<not yet initialised>], latent[<not yet initialised>] 
```
