# Extract model state properties or samples

Extracts model state summary values or generate samples

## Usage

``` r
bru_state(x, ...)

# S3 method for class 'bru_model'
bru_state(x, ..., state = NULL)

# S3 method for class 'inla'
bru_state(
  x,
  property = "mode",
  n = 1,
  seed = 0L,
  num.threads = NULL,
  internal_hyperpar = FALSE,
  ...
)

# S3 method for class 'bru'
bru_state(x, property = "mode", ...)
```

## Arguments

- x:

  Object

- ...:

  Further arguments, passed on to `post.sample.structured()` when
  `property="sample"`.

- state:

  A list of initial state vectors, one for each component. If `NULL`,
  all-zero vectors are returned for each component.

- property:

  Property of the model components to obtain value from. Default:
  "mode". Other options are "mean", "0.025quant", "0.975quant", "sd" and
  "sample", "zeros", "joint_mode", "predictor_sd". In case of "sample"
  you will obtain samples from the posterior. For "zeros", all-zero
  vectors are returned for each component.

- n:

  The number of samples to generate. Default: 1L

- seed:

  Integer seed for random number generation. Default: 0L

- num.threads:

  `num.threads` setting for INLA. Default: NULL, which leaves it up to
  INLA.

- internal_hyperpar:

  logical; If `TRUE`, return hyperparameter properties on the internal
  scale. Currently ignored when `property="sample"`. Default is `FALSE`.

## Value

A list of length `n`

## Examples

``` r
if (bru_safe_inla()) {}
#> NULL
```
