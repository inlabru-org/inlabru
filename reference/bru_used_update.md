# Update used_component information objects

Merge available component labels information with used components
information.

## Usage

``` r
bru_used_update(x, labels)

# S3 method for class 'bru_obs_list'
bru_used_update(x, labels)

# S3 method for class 'bru_obs'
bru_used_update(x, labels, ...)

# S3 method for class 'bru_pred_expr'
bru_used_update(x, labels)

# S3 method for class 'bru_used'
bru_used_update(x, labels, ...)
```

## Arguments

- x:

  Object to be updated

- labels:

  character vector of component labels

## Value

An updated version of `x`. In the `bru_used` information, only
components that are in `labels` are retained. If the ".effect",
".effect.", ".latent", or ".latent." pronoun/container object names are
present in the input `effect` part, all labels will be included in the
output `effect` and `latent` parts, respectively, except for those that
are specifically excluded in an `effect_exclude` part.

## See also

Other bru_used:
[`bru_used()`](https://inlabru-org.github.io/inlabru/reference/bru_used.md),
[`bru_used_vars()`](https://inlabru-org.github.io/inlabru/reference/bru_used_vars.md),
[`new_bru_used()`](https://inlabru-org.github.io/inlabru/reference/new_bru_used.md)
