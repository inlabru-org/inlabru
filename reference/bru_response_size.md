# Response size queries

Extract the number of response values from `bru` and related objects.

## Usage

``` r
bru_response_size(object)

# Default S3 method
bru_response_size(object)

# S3 method for class 'list'
bru_response_size(object)

# S3 method for class 'inla.surv'
bru_response_size(object)

# S3 method for class 'bru_obs'
bru_response_size(object)

# S3 method for class 'bru_obs_list'
bru_response_size(object)

# S3 method for class 'bru_info'
bru_response_size(object)

# S3 method for class 'bru'
bru_response_size(object)
```

## Arguments

- object:

  An object from which to extract response size(s).

## Value

An `integer` vector.

## Methods (by class)

- `bru_response_size(default)`: Extract the number of observations from
  an object supporting [`NROW()`](https://rdrr.io/r/base/nrow.html).

- `bru_response_size(list)`: Extract the number of observations from an
  object supporting `NROW(object[[1]])`.

- `bru_response_size(inla.surv)`: Extract the number of observations
  from an `inla.surv` object.

- `bru_response_size(bru_obs)`: Extract the number of observations from
  a `bru_obs` object.

- `bru_response_size(bru_obs_list)`: Extract the number of observations
  from a `bru_obs_list` object, as a vector with one value per
  observation model.

- `bru_response_size(bru_info)`: Extract the number of observations from
  a `bru_info` object, as a vector with one value per observation model.

- `bru_response_size(bru)`: Extract the number of observations from a
  `bru` object, as a vector with one value per observation model.

## See also

[`like()`](https://inlabru-org.github.io/inlabru/reference/bru_obs.md)

## Examples

``` r
bru_response_size(
  bru_obs(y ~ 1, data = data.frame(y = rnorm(10)), family = "gaussian")
)
#> [1] 10
```
