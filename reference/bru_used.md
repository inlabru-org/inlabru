# List components used in a model

Create or extract information about which components are used by a
model, or its individual observation models. If a non-NULL `labels`
argument is supplied, also calls
[`bru_used_update()`](https://inlabru-org.github.io/inlabru/reference/bru_used_update.md)
on the `bru_used` objects.

## Usage

``` r
bru_used(x = NULL, ...)

# S3 method for class '`NULL`'
bru_used(x = NULL, ...)

# Default S3 method
bru_used(x, ...)

# S3 method for class 'bru'
bru_used(x, ..., join = TRUE)

# S3 method for class 'bru_model'
bru_used(x, ..., join = TRUE)

# S3 method for class 'bru_info'
bru_used(x, ..., join = TRUE)

# S3 method for class 'list'
bru_used(x, ..., join = TRUE)

# S3 method for class 'bru_obs'
bru_used(x, ...)

# S3 method for class 'bru_pred_expr'
bru_used(x, ...)

# S3 method for class 'bru_used'
bru_used(x, labels = NULL, ...)

# S3 method for class 'bru_used'
format(x, ...)

# S3 method for class 'bru_used'
summary(object, ...)

# S3 method for class 'bru_used'
print(x, ...)
```

## Arguments

- x:

  An object that contains information about used components

- ...:

  Parameters passed on to the other methods

- join:

  Whether to join list output into a single object; Default may depend
  on the input object class

- labels:

  character; component labels passed on to
  [`bru_used_update()`](https://inlabru-org.github.io/inlabru/reference/bru_used_update.md)

## Value

A `bru_used` object (a list with elements `effect` and `latent`), or a
list of such objects (for methods with `join = FALSE`)

## Methods (by class)

- `` bru_used(`NULL`) ``: Create a `bru_used` object by calling
  [`new_bru_used()`](https://inlabru-org.github.io/inlabru/reference/new_bru_used.md).

- `bru_used(default)`: Create a `bru_used` object by calling
  [`new_bru_used()`](https://inlabru-org.github.io/inlabru/reference/new_bru_used.md)

- `bru_used(bru)`: Extract the `bru_used` information for the collection
  of observation models used in a `bru` object.

- `bru_used(bru_model)`: Extract the `bru_used` information for the
  collection of observation models used in a `bru_model` object.

- `bru_used(bru_info)`: Extract the `bru_used` information for the
  collection of observation models used in a `bru_info` object.

- `bru_used(list)`: Extract the `bru_used` information for each element
  of a list, and optionally join into a single `bru_used` object.

- `bru_used(bru_obs)`: Extract the `bru_used` information for the
  observation model predictor used in a `bru` observation model
  `bru_obs` object.

- `bru_used(bru_pred_expr)`: Extract the `bru_used` information for an
  observation model predictor `bru_pred_expr` object.

- `bru_used(bru_used)`: Convenience method that takes an existing
  `bru_used` object and calls
  [`bru_used_update()`](https://inlabru-org.github.io/inlabru/reference/bru_used_update.md)
  if `labels` is non-NULL.

## Methods (by generic)

- `format(bru_used)`: Text formatting method for `bru_used` objects.

- `summary(bru_used)`: Summary method for `bru_used` objects.

- `print(bru_used)`: Print method for `bru_used` objects.

## See also

Other bru_used:
[`bru_used_update()`](https://inlabru-org.github.io/inlabru/reference/bru_used_update.md),
[`bru_used_vars()`](https://inlabru-org.github.io/inlabru/reference/bru_used_vars.md),
[`new_bru_used()`](https://inlabru-org.github.io/inlabru/reference/new_bru_used.md)

## Examples

``` r
(used <- new_bru_used(~.))
#> Used effect[<not yet initialised>], latent[]
bru_used(used, labels = c("a", "c"))
#> Used effect[a, c], latent[]
(used <- new_bru_used(~ a + b + c_latent + d_latent))
#> Used effect[a, b], latent[c, d]
bru_used(used, labels = c("a", "c"))
#> Used effect[a], latent[c]
(used <- new_bru_used(~ a_eval(.latent$c)))
#> Used effect[], latent[c, a]
bru_used(used, labels = c("a", "c"))
#> Used effect[], latent[a, c]
```
