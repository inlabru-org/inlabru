# Build an inla data stack from linearisation information

Combine linearisation for multiple likelihoods

## Usage

``` r
bru_make_stack(...)

# S3 method for class 'bru_obs'
bru_make_stack(lhood, lin, idx, ..., family_index = 1L)

# S3 method for class 'bru_obs_list'
bru_make_stack(lhoods, lin, idx, ...)
```

## Arguments

- ...:

  Arguments passed on to other methods

- lhood:

  A `bru_obs` object

- lin:

  Linearisation information

  - For `.bru_obs`, a `bm_taylor` object

  - For `.bru_obs_list`, a list of `bm_taylor` objects

- idx:

  Output from
  [`bru_index.bru_model()`](https://inlabru-org.github.io/inlabru/reference/bru_index.md)

- family_index:

  integer specifying the family sequence index of the observation model

- lhoods:

  A `bru_obs_list` object
