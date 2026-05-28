# Extract a summary property from all results of an inla result

Extract a summary property from all results of an inla result

## Usage

``` r
extract_property(result, property, internal_hyperpar = FALSE)
```

## Arguments

- result:

  an `inla` result object with a `bru_info` element

- property:

  character; "mean", "sd", "mode", or some other column identifier for
  inla result `$summary.fixed`, `$summary.random$label`, and
  `$summary.hyperpar`, or "joint_mode". For "joint_mode", the joint
  latent mode is extracted, and the joint hyperparameter mode, in the
  internal scale. For "predictor_sd" the posterior standard deviations
  of the linear predictor are returned.

- internal_hyperpar:

  logical; if `TRUE`, use internal scale for hyperparamter properties.
  Default is `FALSE`, except when `property` is "joint_mode" which
  forces `internal_hyperpar=TRUE`.

## Value

named list for each estimated fixed effect coefficient, random effect
vector, and hyperparameter. The hyperparameter names are standardised
with
[`bru_standardise_names()`](https://inlabru-org.github.io/inlabru/reference/bru_standardise_names.md)
