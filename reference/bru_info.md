# Methods for bru_info objects

The `bru_info` class is used to store metadata about `bru` models.

## Usage

``` r
bru_info(...)

# S3 method for class 'character'
bru_info(method, ..., inlabru_version = NULL, INLA_version = NULL)

# S3 method for class 'bru'
bru_info(object, ...)

as_bru_info(object, ...)

# S3 method for class 'bru_info'
as_bru_info(object, ...)

# S3 method for class 'bru'
as_bru_info(object, ...)

as_bru_model(object, ...)

# S3 method for class 'bru_model'
as_bru_model(object, ...)

# S3 method for class 'bru_info'
as_bru_model(object, ...)

# S3 method for class 'bru'
as_bru_model(object, ...)

# S3 method for class 'bru_info'
summary(object, verbose = TRUE, ...)

# S3 method for class 'summary_bru_info'
print(x, ...)

# S3 method for class 'bru_info'
print(x, ...)
```

## Arguments

- ...:

  Additional arguments to be stored in the `bru_info` object. For
  `summary` and `print` methods, arguments passed on to submethods.

- method:

  character; The type of estimation method used

- inlabru_version:

  character; inlabru package version. Default: NULL, for automatically
  detecting the version

- INLA_version:

  character; INLA package version. Default: NULL, for automatically
  detecting the version

- object:

  A `bru_info` object

- verbose:

  logical; If `TRUE`, include more details of the component definitions.
  If `FALSE`, only show basic component definition information. Default:
  `FALSE`

- x:

  An object to be printed

## Value

A `bru_info` object

## Methods (by class)

- `bru_info(character)`: Create a `bru_info` object

- `bru_info(bru)`: Extract the `bru_info` object from an estimated
  [`bru()`](https://inlabru-org.github.io/inlabru/reference/bru.md)
  result object. The default print method shows information about model
  components and observation models.

## Methods (by generic)

- `as_bru_info(bru_info)`: Extract a `bru_info` object.

- `as_bru_model(bru_info)`: Extract a `bru_model` object.

## Functions

- `as_bru_info()`: Extract the `bru_info` object from an estimated
  [`bru()`](https://inlabru-org.github.io/inlabru/reference/bru.md)
  result object. The default print method shows information about model
  components and observation models.

- `as_bru_info(bru)`: Extract the `bru_info` object from an estimated
  [`bru()`](https://inlabru-org.github.io/inlabru/reference/bru.md)
  result object.

- `as_bru_model()`: Extract the `bru_model` object from an estimated
  [`bru()`](https://inlabru-org.github.io/inlabru/reference/bru.md)
  result object.

- `as_bru_model(bru_model)`: Extract a `bru_model` object.

- `as_bru_model(bru)`: Extract the `bru_model` object from an estimated
  [`bru()`](https://inlabru-org.github.io/inlabru/reference/bru.md)
  result object.
