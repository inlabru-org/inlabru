# Obtain bru component evaluation functions

Internal generic function and methods for generating evaluation
functions for `bru_comp` objects, that can be used in
[`generate.bru()`](https://inlabru-org.github.io/inlabru/reference/generate.md)
and
[`predict.bru()`](https://inlabru-org.github.io/inlabru/reference/predict.md)
expressions.

## Usage

``` r
bru_eval_fun(x, ...)

# S3 method for class 'bru_comp'
bru_eval_fun(x, ...)

# S3 method for class 'bru_comp_list'
bru_eval_fun(x, ...)

# S3 method for class 'bru_model'
bru_eval_fun(x, ...)
```

## Value

A component evaluation function or a list of functions

## Methods (by class)

- `bru_eval_fun(bru_comp)`: Returns function that takes the arguments
  `main`, `group`, `replicate`, `weights`, and `.state` and can be
  evaluated in the context of a data mask from
  [`bru_data_mask()`](https://inlabru-org.github.io/inlabru/reference/bru_data_mask.md).

- `bru_eval_fun(bru_comp_list)`: Returns a list of component evaluation
  functions named `<label>_eval` for each component with label
  `<label>`.

- `bru_eval_fun(bru_model)`: Returns a list of component evaluation
  functions named `<label>_eval` for each model component with label
  `<label>`.
