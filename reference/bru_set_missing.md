# Set missing values in observation models

Set all or parts of the observation model response data to `NA`, for
example for use in cross validation (with
[`bru_rerun()`](https://inlabru-org.github.io/inlabru/reference/bru.md))
or prior sampling (with
[`bru_rerun()`](https://inlabru-org.github.io/inlabru/reference/bru.md)
and
[`generate()`](https://inlabru-org.github.io/inlabru/reference/generate.md),
but see "Prior sampling caveats" below).

## Usage

``` r
bru_set_missing(object, keep = FALSE, ...)

bru_set_missing(x, ...) <- value

# S3 method for class 'bru'
bru_set_missing(object, keep = FALSE, ...)

# S3 method for class 'bru_info'
bru_set_missing(object, keep = FALSE, ...)

# S3 method for class 'bru_obs_list'
bru_set_missing(object, keep = FALSE, ...)

# S3 method for class 'bru_obs'
bru_set_missing(object, keep = FALSE, ...)

# Default S3 method
bru_set_missing(object, keep = FALSE, ...)

# S3 method for class 'data.frame'
bru_set_missing(object, keep = FALSE, ...)

# S3 method for class 'inla.surv'
bru_set_missing(object, keep = FALSE, ...)
```

## Arguments

- object:

  A `bru`, `bru_obs` or `bru_obs_list` object

- keep:

  For `bru_obs`, a single logical or an integer vector; If `TRUE`, keep
  all the response data, if `FALSE` (default), set all of it to `NA`. An
  integer vector determines which elements to keep (for positive values)
  or to set as missing (negative values).

  For `bru` and `bru_obs_list`, a logical scalar or vector, or a list,
  see Details.

- ...:

  Additional arguments passed on to the `bru_obs` and individual data
  class methods. Currently unused.

- x:

  Object on which to apply `bru_set_missing()`

- value:

  Value to be passed as `keep` to `bru_set_missing()`

## Details

For `bru` and `bru_obs_list`,

- `keep` must be either a single logical, which is expanded to a list,

- a logical vector, which is converted to a list,

- an unnamed list of the same length as the number of observation
  models, with elements compatible with the `bru_obs` method, or

- a named list with elements compatible with the `bru_obs` method, and
  only the named `bru_obs` models are acted upon, i.e. the elements not
  present in the list are treated as `keep = TRUE`.

E.g.: `keep = list(b = FALSE)` sets all observations in model `b` to
missing, and does not change model `a`.

E.g.: `keep = list(a = 1:4, b = -(3:5))` keeps only observations `1:4`
of model `a`, marking the rest as missing, and sets observations `3:5`
of model `b` to missing.

## Methods (by class)

- `bru_set_missing(default)`: From `> 2.13.0`, set missing values in any
  object supporting [`base::is.na<-()`](https://rdrr.io/r/base/NA.html)
  for positive and negative indices.

- `bru_set_missing(data.frame)`: From `> 2.13.0`, handles `data.frame`,
  tibbles, including `inla.mdata`.

- `bru_set_missing(inla.surv)`: From `> 2.13.0`, handles `inla.surv`.

## Functions

- `bru_set_missing(x, ...) <- value`: Setter method for
  `bru_set_missing()`

## Prior sampling caveats

Note that prior sampling requires special care for hyperparameters, as
the prior modes are not typically useful; in the future, we plan to have
a dedicated method that samples from the hyperparameters, and then uses

    bru_rerun(
      bru_set_missing(...),
      options = list(
        control.mode = list(
          theta = theta_sample,
          fixed = TRUE)
        )
      )

for each sample.

## Examples

``` r
obs <- c(
  A = bru_obs(y_A ~ ., data = data.frame(y_A = 1:6)),
  B = bru_obs(y_B ~ ., data = data.frame(y_B = 11:15))
)
bru_response_size(obs)
#> A B 
#> 6 5 
lapply(
  bru_set_missing(obs, keep = FALSE),
  function(x) {
    x[["response_data"]][[x[["response"]]]]
  }
)
#> $A
#> [1] NA NA NA NA NA NA
#> 
#> $B
#> [1] NA NA NA NA NA
#> 
lapply(
  bru_set_missing(obs, keep = list(B = FALSE)),
  function(x) {
    x[["response_data"]][[x[["response"]]]]
  }
)
#> $A
#> [1] 1 2 3 4 5 6
#> 
#> $B
#> [1] NA NA NA NA NA
#> 
lapply(
  bru_set_missing(obs, keep = list(1:4, -(3:5))),
  function(x) {
    x[["response_data"]][[x[["response"]]]]
  }
)
#> $A
#> [1]  1  2  3  4 NA NA
#> 
#> $B
#> [1] 11 12 NA NA NA
#> 

(obs <- INLA::inla.mdata(y = 1:4, X = matrix(1:8, 4, 2)))
#> inla.cols =  2 1 2 
#>   Y1 X1 X2
#> 1  1  1  5
#> 2  2  2  6
#> 3  3  3  7
#> 4  4  4  8
bru_set_missing(obs, keep = c(1, 4))
#> inla.cols =  2 1 2 
#>   Y1 X1 X2
#> 1  1  1  5
#> 2 NA NA NA
#> 3 NA NA NA
#> 4  4  4  8
bru_set_missing(obs) <- -(1:2)
obs
#> inla.cols =  2 1 2 
#>   Y1 X1 X2
#> 1 NA NA NA
#> 2 NA NA NA
#> 3  3  3  7
#> 4  4  4  8
(obs <- INLA::inla.surv(time = 1:4, event = c(1, 0, 1, 0)))
#> [1] 1  2+ 3  4+
bru_set_missing(obs, keep = c(1, 4))
#> [1] 1   NA+ NA  4+ 

(obs <- INLA::inla.surv(
  time = 1:4,
  event = c(1, 0, 1, 0),
  cure = matrix(1:8, 4, 2)
))
#> [1] 1  2+ 3  4+
bru_set_missing(obs, keep = c(1, 4))
#> [1] 1   NA+ NA  4+ 

(obs <- INLA::inla.surv(
  time = 1:4,
  event = c(1, 0, 1, 0),
  subject = c(1, 1, 2, 1)
))
#> [1] 1  0+ 3  0+
bru_set_missing(obs, keep = c(1, 4))
#> [1] 1   NA+ NA  0+ 
```
