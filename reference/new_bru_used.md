# Store information about used components

Create a `bru_used` object from effect name character vectors. Create
information about which components are used by a model, or its
individual observation models. If a non-NULL `labels` argument is
supplied, also calls
[`bru_used_update()`](https://inlabru-org.github.io/inlabru/reference/bru_used_update.md)
on the `bru_used` object.

## Usage

``` r
new_bru_used(
  x = NULL,
  ...,
  effect = NULL,
  effect_exclude = NULL,
  latent = NULL,
  labels = NULL
)

# S3 method for class '`NULL`'
new_bru_used(
  x = NULL,
  ...,
  effect = NULL,
  effect_exclude = NULL,
  latent = NULL,
  labels = NULL
)

# Default S3 method
new_bru_used(
  x,
  ...,
  effect = NULL,
  effect_exclude = NULL,
  latent = NULL,
  labels = NULL
)

# S3 method for class 'character'
new_bru_used(
  x,
  ...,
  effect = NULL,
  effect_exclude = NULL,
  latent = NULL,
  labels = NULL
)
```

## Arguments

- x:

  `NULL`, or an object representing an expression

- ...:

  Parameters passed on to the other methods

- effect:

  character; components used as effects. When `NULL`, auto-detect
  components to include all components in a predictor expression.

- effect_exclude:

  character; components to specifically exclude from effect evaluation.
  When `NULL`, do not specifically exclude any components.

- latent:

  character; components used as `_latent` or `_eval()`. When `NULL`,
  auto-detect components.

- labels:

  character; component labels passed on to
  [`bru_used_update()`](https://inlabru-org.github.io/inlabru/reference/bru_used_update.md)

## Value

A `bru_used` object (a list with elements `effect` and `latent`)

## Details

The arguments `effect`, `effect_exclude`, and `latent` control what
components and effects are available for use in predictor expressions.

- `effect`:

  Character vector of component labels that are used as effects by the
  predictor expression; If `NULL` (default), the names are extracted
  from the formula.

- `exclude`:

  Character vector of component labels to be excluded from the effect
  list even if they have been auto-detected as being necessary. Default
  is `NULL`; do not remove any components from the inclusion list.

- `include_latent`:

  Character vector. Specifies which latent state variables need to be
  directly available to the predictor expression, with a `_latent`
  suffix. This also makes evaluator functions with suffix `_eval`
  available, taking parameters `main`, `group`, and `replicate`, taking
  values for where to evaluate the component effect that are different
  than those defined in the component definition itself (see
  [`bru_comp_eval()`](https://inlabru-org.github.io/inlabru/reference/bru_comp_eval.md)).
  If `NULL`, the use of `_latent` and `_eval` in the predictor
  expression is detected automatically.

## Methods (by class)

- `` new_bru_used(`NULL`) ``: Create a `bru_used` object from effect
  name character vectors.

- `new_bru_used(default)`: Create a `bru_used` object from an expression
  object supported by
  [`bru_used_vars()`](https://inlabru-org.github.io/inlabru/reference/bru_used_vars.md).

- `new_bru_used(character)`: Create a `bru_used` object from a string.

## See also

Other bru_used:
[`bru_used()`](https://inlabru-org.github.io/inlabru/reference/bru_used.md),
[`bru_used_update()`](https://inlabru-org.github.io/inlabru/reference/bru_used_update.md),
[`bru_used_vars()`](https://inlabru-org.github.io/inlabru/reference/bru_used_vars.md)

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
```
