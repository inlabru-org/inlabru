# Extract predictor or component index information

**\[experimental\]** Extract the index vector for a
[`bru_obs()`](https://inlabru-org.github.io/inlabru/reference/bru_obs.md)
predictor, or the whole or a subset of a full
[`bru()`](https://inlabru-org.github.io/inlabru/reference/bru.md)
predictor.

## Usage

``` r
bru_index(object, ...)

# S3 method for class 'bru_obs'
bru_index(object, what = NULL, ...)

# S3 method for class 'bru_obs_list'
bru_index(object, tag = NULL, what = NULL, ...)

# S3 method for class 'bru'
bru_index(object, tag = NULL, what = NULL, ...)

# S3 method for class 'bru_info'
bru_index(object, tag = NULL, what = NULL, ...)

# S3 method for class 'bru_model'
bru_index(object, used, ...)

# S3 method for class 'bru_comp'
bru_index(object, inla_f, ...)

# S3 method for class 'bru_comp_list'
bru_index(object, inla_f, ...)

index_eval(...)
```

## Arguments

- object:

  A
  [component](https://inlabru-org.github.io/inlabru/reference/bru_comp.md).

- ...:

  Arguments passed on to sub-methods.

- what:

  `character` or `NULL`; One of `NULL`, "all", "observed", and
  "missing". If `NULL` (default) or "all", gives the index vector for
  the full sub-model predictor. If "observed", gives the index vector
  for the observed part (response is not `NA`). If "missing", gives the
  index vector for the missing part (response is `NA`) of the model.

- tag:

  `character` or `integer`; Either a character vector identifying the
  tags of one or more of the
  [`bru_obs()`](https://inlabru-org.github.io/inlabru/reference/bru_obs.md)
  observation models, or an integer vector identifying models by their
  [`bru()`](https://inlabru-org.github.io/inlabru/reference/bru.md)
  specification order. If `NULL` (default) computes indices for all
  sub-models.

- used:

  A
  [`bru_used()`](https://inlabru-org.github.io/inlabru/reference/bru_used.md)
  object

- inla_f:

  logical; when `TRUE`, must result in values compatible with
  `INLA::f(...)` an specification and corresponding
  `INLA::inla.stack(...)` constructions.

## Value

- `bru_index(bru_obs)`: An `integer` vector.

&nbsp;

- `bru_index(bru_obs_list)`: An `integer` vector.

&nbsp;

- `bru_index(bru)`: An `integer` vector.

&nbsp;

- `bru_index(bru_info)`: An `integer` vector.

&nbsp;

- `bru_index(bru_model)`: A named list of `idx_full` and `idx_inla`,
  named list of indices, and `inla_subset`, and `inla_subset`, a named
  list of logical subset specifications for extracting the
  [`INLA::f()`](https://rdrr.io/pkg/INLA/man/f.html) compatible index
  subsets.

&nbsp;

- `bru_index(bru_comp)`: A list of indices into the latent variables
  compatible with the component mapper.

&nbsp;

- `bru_index(bru_comp_list)`: A list of list of indices into the latent
  variables compatible with each component mapper.

## Methods (by class)

- `bru_index(bru_obs)`: Extract the index vector for the predictor
  vector for a
  [`bru_obs()`](https://inlabru-org.github.io/inlabru/reference/bru_obs.md)
  sub-model. The indices are relative to the sub-model, and need to be
  appropriately offset to be used in the full model predictor.

- `bru_index(bru_obs_list)`: Extract the index vector for "APredictor"
  for one or more specified observation
  [`bru_obs()`](https://inlabru-org.github.io/inlabru/reference/bru_obs.md)
  sub-models. Accepts any combination of `tag` and `what`.

- `bru_index(bru)`: Extract the index vector for "APredictor" for one or
  more specified observation
  [`bru_obs()`](https://inlabru-org.github.io/inlabru/reference/bru_obs.md)
  sub-models. Accepts any combination of `tag` and `what`.

- `bru_index(bru_info)`: Extract the index vector for "APredictor" for
  one or more specified observation
  [`bru_obs()`](https://inlabru-org.github.io/inlabru/reference/bru_obs.md)
  sub-models. Accepts any combination of `tag` and `what`.

- `bru_index(bru_model)`: Compute all index values for a
  [`bru_model()`](https://inlabru-org.github.io/inlabru/reference/bru_model.md)
  object. Computes the index value matrices for included components
  according to the `used` argument.

## Functions

- `index_eval()`: **\[deprecated\]** since "2.12.0.9023". Use
  `bru_index()` instead.

## Author

Fabian E. Bachl <bachlfab@gmail.com>, Finn Lindgren
<finn.lindgren@gmail.com>

## Examples

``` r
fit <- bru(
  ~ 0 + x,
  bru_obs(
    y ~ .,
    data = data.frame(x = 1:3, y = 1:3 + rnorm(3)),
    tag = "A"
  ),
  bru_obs(
    y ~ .,
    data = data.frame(x = 1:4, y = c(NA, NA, 3:4) + rnorm(4)),
    tag = "B"
  ),
  options = list(bru_run = FALSE) # We only need the model structure
)
bru_index(fit)
#> [1] 1 2 3 4 5 6 7
bru_index(fit, "A")
#> [1] 1 2 3
bru_index(fit, "B")
#> [1] 4 5 6 7
bru_index(fit, c("B", "A"))
#> [1] 4 5 6 7 1 2 3
bru_index(fit, what = "missing")
#> [1] 4 5
```
