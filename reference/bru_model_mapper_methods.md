# Mapper methods for model objects

Methods for the
[`ibm_as_taylor()`](https://inlabru-org.github.io/inlabru/reference/ibm_as_taylor.md)
and
[`ibm_simplify()`](https://inlabru-org.github.io/inlabru/reference/ibm_simplify.md)
methods for
[bru](https://inlabru-org.github.io/inlabru/reference/bru.md) model
objects and related classes.

## Usage

``` r
# S3 method for class 'bru_model'
ibm_as_taylor(
  mapper,
  input,
  state = NULL,
  ...,
  options = NULL,
  comp_mappers = NULL,
  eval_fun = NULL
)

# S3 method for class 'bru_comp_list'
ibm_as_taylor(mapper, input, state = NULL, ...)

# S3 method for class 'bru_model'
ibm_simplify(mapper, input = NULL, state = NULL, ...)

# S3 method for class 'bru_comp'
ibm_simplify(mapper, input = NULL, state = NULL, ...)

# S3 method for class 'bru_comp_list'
ibm_simplify(mapper, input = NULL, state = NULL, ...)

# S3 method for class 'bm_list'
ibm_as_taylor(mapper, input, state = NULL, ...)

# S3 method for class 'bm_list'
ibm_simplify(mapper, input = NULL, state = NULL, ...)

# S3 method for class 'bm_list'
ibm_eval2(mapper, input, state = NULL, ...)

# S3 method for class 'bm_list'
ibm_eval(mapper, input, state = NULL, ...)

# S3 method for class 'bm_list'
ibm_jacobian(mapper, input, state = NULL, ...)

# S3 method for class 'bru_model'
ibm_eval2(
  mapper,
  input,
  state = NULL,
  ...,
  options = NULL,
  comp_mappers = NULL,
  eval_fun = NULL
)

# S3 method for class 'bru_comp_list'
ibm_eval(mapper, input, state = NULL, ..., comp_mappers = NULL)

# S3 method for class 'bru_comp'
ibm_eval(mapper, input, state = NULL, ...)

# S3 method for class 'bru_model'
ibm_eval(
  mapper,
  input,
  state = NULL,
  ...,
  options = NULL,
  comp_mappers = NULL,
  eval_fun = NULL
)

# S3 method for class 'bru_model'
ibm_jacobian(
  mapper,
  input,
  state = NULL,
  ...,
  options = NULL,
  comp_mappers = NULL,
  eval_fun = NULL
)
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

- options:

  A
  [bru_options](https://inlabru-org.github.io/inlabru/reference/bru_options.md)
  options object or a list of options passed on to
  [`bru_options()`](https://inlabru-org.github.io/inlabru/reference/bru_options.md)

- comp_mappers:

  A
  [`bm_list()`](https://inlabru-org.github.io/inlabru/reference/bm_list.md)
  of mappers, either the original component mappers, or simplified
  mappers.

- eval_fun:

  A list of functions, typically from
  [`bru_eval_fun()`](https://inlabru-org.github.io/inlabru/reference/bru_eval_fun.md).

## Functions

- `ibm_as_taylor(bru_model)`: Returns a list (one element per
  observation model) of
  [bm_list](https://inlabru-org.github.io/inlabru/reference/bm_list.md)
  objects, each with one
  [bm_taylor](https://inlabru-org.github.io/inlabru/reference/bm_taylor.md)
  entry for each included component. (autodiff == "pandemic")

  If `autodiff` is not "pandemic", returns a list of
  [bm_taylor](https://inlabru-org.github.io/inlabru/reference/bm_taylor.md)
  objects, one for each observation model, with the offset and jacobians
  evaluated for the predictor of the observation model, and the
  component mappers passed on as `comp_mappers` for the evaluation of
  the jacobians.

- `ibm_simplify(bru_model)`: Returns a list (one element per observation
  model) of
  [bm_list](https://inlabru-org.github.io/inlabru/reference/bm_list.md)
  objects, each with one
  [bru_mapper](https://inlabru-org.github.io/inlabru/reference/bru_mapper.md)
  entry for each included component.
