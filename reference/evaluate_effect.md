# Evaluate a component effect

Calculate latent component effects given some data and the state of the
component's internal random variables.

## Usage

``` r
evaluate_effect_single_state(...)

# S3 method for class 'bru_mapper'
evaluate_effect_single_state(component, input, state, ..., label = NULL)

# S3 method for class 'bm_list'
evaluate_effect_single_state(components, input, state, ...)

# S3 method for class 'bru_comp_list'
evaluate_effect_single_state(components, input, state, ...)
```

## Arguments

- ...:

  Optional additional parameters, e.g. `inla_f`. Normally unused.

- component:

  A
  [bru_mapper](https://inlabru-org.github.io/inlabru/reference/bru_mapper.md),
  [bru_comp](https://inlabru-org.github.io/inlabru/reference/bru_comp.md),
  or
  [bm_list](https://inlabru-org.github.io/inlabru/reference/bm_list.md).

- input:

  Pre-evaluated component input

- state:

  Specification of one latent variable state:

  - `evaluate_effect_single_state.bru_mapper`: A vector of the latent
    component state.

  - `evaluate_effect_single_state.*_list`: list of named state vectors.

- label:

  Option label used for any warning messages, specifying the affected
  component.

## Value

- `evaluate_effect_single_state.bm_list`: A list of evaluated component
  effect values

## Author

Fabian E. Bachl <bachlfab@gmail.com> and Finn Lindgren
<finn.lindgren@gmail.com>
