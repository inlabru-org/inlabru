# Mapper methods for model objects

Methods for the
[`ibm_linear()`](https://inlabru-org.github.io/inlabru/reference/ibm_linear.md)
and
[`ibm_simplify()`](https://inlabru-org.github.io/inlabru/reference/ibm_simplify.md)
methods for
[bru](https://inlabru-org.github.io/inlabru/reference/bru.md) model
objects and related classes.

## Usage

``` r
# S3 method for class 'bru_model'
ibm_linear(mapper, input, state = NULL, ...)

# S3 method for class 'bru_comp_list'
ibm_linear(mapper, input, state = NULL, ...)

# S3 method for class 'bru_model'
ibm_simplify(mapper, input = NULL, state = NULL, ...)

# S3 method for class 'bru_comp'
ibm_simplify(mapper, input = NULL, state = NULL, ...)

# S3 method for class 'bru_comp_list'
ibm_simplify(mapper, input = NULL, state = NULL, ...)

# S3 method for class 'bm_list'
ibm_linear(mapper, input, state = NULL, ...)

# S3 method for class 'bm_list'
ibm_simplify(mapper, input = NULL, state = NULL, ...)
```

## Arguments

- mapper:

  A mapper S3 object, inheriting from `bru_mapper`.

- input:

  Data input for the mapper.

- state:

  A vector of latent state values for the mapping, of length
  `ibm_n(mapper, inla_f = FALSE)`

- ...:

  Arguments passed on to other methods

## Functions

- `ibm_linear(bru_model)`: Returns a list (one element per observation
  model) of
  [bm_list](https://inlabru-org.github.io/inlabru/reference/bm_list.md)
  objects, each with one
  [bm_taylor](https://inlabru-org.github.io/inlabru/reference/bm_taylor.md)
  entry for each included component.

- `ibm_simplify(bru_model)`: Returns a list (one element per observation
  model) of
  [bm_list](https://inlabru-org.github.io/inlabru/reference/bm_list.md)
  objects, each with one
  [bru_mapper](https://inlabru-org.github.io/inlabru/reference/bru_mapper.md)
  entry for each included component.
