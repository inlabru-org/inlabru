# Extract mapper information from INLA model component objects

The component definitions will automatically attempt to extract mapper
information from any model object by calling the generic
`bru_get_mapper`. Any class method implementation should return a
[bru_mapper](https://inlabru-org.github.io/inlabru/reference/bru_mapper.md)
object suitable for the given latent model.

## Usage

``` r
bru_get_mapper(model, ...)

# S3 method for class 'inla.spde'
bru_get_mapper(model, ...)

# S3 method for class 'inla.rgeneric'
bru_get_mapper(model, ...)

# S3 method for class 'inla.cgeneric'
bru_get_mapper(model, ...)

bru_get_mapper_safely(model, ...)

# S3 method for class 'inla_model_reparam'
bru_get_mapper(model, ...)
```

## Arguments

- model:

  A model component object

- ...:

  Arguments passed on to other methods

## Value

A
[bru_mapper](https://inlabru-org.github.io/inlabru/reference/bru_mapper.md)
object defined by the model component

## Details

Before implementing your own `bru_get_mapper` method, check if there is
already a general method available that handles your model class, such
as `bru_get_mapper.inla.spde()`, `bru_get_mapper.inla.rgeneric()`, and
`bru_get_mapper.inla.cgeneric()`.

## Methods (by class)

- `bru_get_mapper(inla.spde)`: Extract an indexed mapper for the
  `model$mesh` object contained in the model object, which is assumed to
  be of a class supporting relevant `fmesher` methods.

- `bru_get_mapper(inla.rgeneric)`: Returns the mapper given by a
  pre-computed mapper, or an index mapper mapping the size of the model
  graph.

  The easiest method to define a mapper for an `inla.rgeneric` model is
  to store the mapper in the object. Alternatively, define your model
  using a subclass and define a `bru_get_mapper.<subclass>` method that
  should return the corresponding `bru_mapper` object.

  The order of precedence for the mapper construction when calling
  `bru_get_mapper(model)` has the following precedence:

  1.  `bru_get_mapper.<subclass>`, if `model` has a subclass, otherwise

  2.  `model[["mapper"]]` if that is `NULL`, and otherwise

  3.  [`bm_index()`](https://inlabru-org.github.io/inlabru/reference/bm_index.md)
      using the size of the graph returned by `model[["f"]][["n"]]`

- `bru_get_mapper(inla.cgeneric)`: Works the same as the method of
  `inla.rgeneric`

- `bru_get_mapper(inla_model_reparam)`: Reparameterised inla model
  mapper, see
  [`inla.spde2.pcmatern_B()`](https://inlabru-org.github.io/inlabru/reference/pcmatern_B.md)
  **\[experimental\]**

## Functions

- `bru_get_mapper_safely()`: Tries to call the `bru_get_mapper`, and
  returns `NULL` if it fails (e.g. due to no available class method). If
  the call succeeds and returns non-`NULL`, it checks that the object
  inherits from the `bru_mapper` class, and gives an error if it does
  not.

## See also

[bru_mapper](https://inlabru-org.github.io/inlabru/reference/bru_mapper.md)
for mapper constructor methods, and the individual mappers for specific
implementation details.

Other mappers:
[`bm_aggregate()`](https://inlabru-org.github.io/inlabru/reference/bm_aggregate.md),
[`bm_collect()`](https://inlabru-org.github.io/inlabru/reference/bm_collect.md),
[`bm_const()`](https://inlabru-org.github.io/inlabru/reference/bm_const.md),
[`bm_expr()`](https://inlabru-org.github.io/inlabru/reference/bm_expr.md),
[`bm_factor()`](https://inlabru-org.github.io/inlabru/reference/bm_factor.md),
[`bm_fm_mesh_1d`](https://inlabru-org.github.io/inlabru/reference/bm_fm_mesh_1d.md),
[`bm_fmesher()`](https://inlabru-org.github.io/inlabru/reference/bm_fmesher.md),
[`bm_harmonics()`](https://inlabru-org.github.io/inlabru/reference/bm_harmonics.md),
[`bm_index()`](https://inlabru-org.github.io/inlabru/reference/bm_index.md),
[`bm_linear()`](https://inlabru-org.github.io/inlabru/reference/bm_linear.md),
[`bm_logitaverage()`](https://inlabru-org.github.io/inlabru/reference/bm_logitaverage.md),
[`bm_logsumexp()`](https://inlabru-org.github.io/inlabru/reference/bm_logsumexp.md),
[`bm_marginal()`](https://inlabru-org.github.io/inlabru/reference/bm_marginal.md),
[`bm_matrix()`](https://inlabru-org.github.io/inlabru/reference/bm_matrix.md),
[`bm_multi()`](https://inlabru-org.github.io/inlabru/reference/bm_multi.md),
[`bm_pipe()`](https://inlabru-org.github.io/inlabru/reference/bm_pipe.md),
[`bm_reparam()`](https://inlabru-org.github.io/inlabru/reference/bm_reparam.md),
[`bm_repeat()`](https://inlabru-org.github.io/inlabru/reference/bm_repeat.md),
[`bm_scale()`](https://inlabru-org.github.io/inlabru/reference/bm_scale.md),
[`bm_shift()`](https://inlabru-org.github.io/inlabru/reference/bm_shift.md),
[`bm_sum()`](https://inlabru-org.github.io/inlabru/reference/bm_sum.md),
[`bm_taylor()`](https://inlabru-org.github.io/inlabru/reference/bm_taylor.md),
[`bru_mapper()`](https://inlabru-org.github.io/inlabru/reference/bru_mapper.md)

## Examples

``` r
if (bru_safe_inla()) {
  mesh <- fmesher::fm_rcdt_2d_inla(globe = 2)
  spde <- INLA::inla.spde2.pcmatern(mesh,
    prior.range = c(1, 0.5),
    prior.sigma = c(1, 0.5)
  )
  mapper <- bru_get_mapper(spde)
  ibm_n(mapper)
}
#> [1] 42
```
