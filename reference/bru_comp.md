# Latent model component construction

Similar to [`glm()`](https://rdrr.io/r/stats/glm.html), `gam()` and
[`inla()`](https://rdrr.io/pkg/INLA/man/inla.html),
[`bru()`](https://inlabru-org.github.io/inlabru/reference/bru.md) models
can be constructed via a formula-like syntax, where each latent effect
is specified. However, in addition to the parts of the syntax compatible
with [`INLA::inla`](https://rdrr.io/pkg/INLA/man/inla.html), `bru`
components offer additional functionality which facilitates modelling,
and the predictor expression can be specified separately, allowing more
complex and non-linear predictors to be defined. The formula syntax is
just a way to allow all model components to be defined in a single line
of code, but the definitions can optionally be split up into separate
component definitions. See Details for more information.

The `bru_comp` methods all rely on the `bru_comp.character()` method,
that defines a model component with a given label/name. The user usually
doesn't need to call these methods directly, but can instead supply a
formula expression that can be interpreted by the
[`bru_comp_list.formula()`](https://inlabru-org.github.io/inlabru/reference/bru_comp_list.md)
method, called inside
[`bru()`](https://inlabru-org.github.io/inlabru/reference/bru.md).

## Usage

``` r
bru_comp(...)

bru_component(...)

# S3 method for class 'character'
bru_comp(
  object,
  main = NULL,
  weights = NULL,
  ...,
  model = NULL,
  mapper = NULL,
  main_layer = NULL,
  main_selector = NULL,
  n = NULL,
  values = NULL,
  season.length = NULL,
  nrow = NULL,
  ncol = NULL,
  copy = NULL,
  weights_layer = NULL,
  weights_selector = NULL,
  group = 1L,
  group_mapper = NULL,
  group_layer = NULL,
  group_selector = NULL,
  ngroup = NULL,
  control.group = NULL,
  replicate = 1L,
  replicate_mapper = NULL,
  replicate_layer = NULL,
  replicate_selector = NULL,
  nrep = NULL,
  marginal = NULL,
  A.msk = deprecated(),
  .envir = parent.frame(),
  envir_extra = NULL
)
```

## Arguments

- ...:

  Parameters passed on to other methods

- object:

  A character label for the component

- main:

  `main` takes an R expression that evaluates to where the latent
  variables should be evaluated (coordinates, indices, continuous scalar
  (for rw2 etc)). Arguments starting with weights, group, replicate
  behave similarly to main, but for the corresponding features of
  [`INLA::f()`](https://rdrr.io/pkg/INLA/man/f.html). `main` supports
  `tidy` evaluation expressions, with `.data` and `.env` pronouns.

- weights, weights_layer, weights_selector:

  Optional specification of effect scaling weights. Same syntax as for
  `main`.

- model:

  Either one of "const" (same as "offset"), "factor_full",
  "factor_contrast", "linear", "fixed", or a model name or object
  accepted by INLA's `f` function. If set to NULL, then "linear" is used
  for vector inputs, and "fixed" for matrix input (converted internally
  to an iid model with fixed precision)

- mapper:

  Information about how to do the mapping from the values evaluated in
  `main`, and to the latent variables. Auto-detects spde model objects
  in model and extracts the mesh object to use as the mapper, and
  auto-generates mappers for indexed models. (Default: NULL, for
  auto-determination)

- main_layer, main_selector:

  The `_layer` input should evaluate to a numeric index or character
  name or vector of which layer/variable to extract from a covariate
  data object given in `main`. (Default: NULL if `_selector` is given.
  Otherwise the effect component name, if it exists in the covariate
  object, and otherwise the first column of the covariate data frame)

  The `_selector` value should be a character name of a variable whose
  contents determines which layer to extract from a covariate for each
  data point. (Default: NULL)

- n:

  The number of latent variables in the model. Should be auto-detected
  for most or all models. Default: NULL, for auto-detection. Models with
  matrix input for `Cmatrix` or `graph` will use the matrix size. An
  error is given if it realises it can't figure it out by itself.

- values:

  Specifies for what covariate/index values INLA should build the latent
  model. Normally generated internally based on the mapping details.
  (Default: NULL, for auto-determination)

- season.length:

  Passed on to [`INLA::f()`](https://rdrr.io/pkg/INLA/man/f.html) for
  model `"seasonal"` (TODO: check if this parameter is still fully
  handled)

- nrow, ncol:

  Number of rows and columns for `model` types "rw2d", "rw2diid", and
  "matern2d". Default is `NULL`.

- copy:

  character; label of other component that this component should be a
  copy of. If the `fixed = FALSE`, a scaling constant is estimated, via
  a hyperparameter. If `fixed = TRUE`, the component scaling is fixed,
  by default to 1; for fixed scaling, it's more efficient to express the
  scaling in the predictor expression instead of making a copy
  component.

- group, group_mapper, group_layer, group_selector, ngroup:

  Optional specification of kronecker/group model indexing.

- control.group:

  `list` of kronecker/group model parameters, currently passed directly
  on to [`INLA::f`](https://rdrr.io/pkg/INLA/man/f.html)

- replicate, replicate_mapper, replicate_layer, replicate_selector,
  nrep:

  Optional specification of indices for an independent replication
  model. Same syntax as for `main`

- marginal:

  May specify a
  [`bm_marginal()`](https://inlabru-org.github.io/inlabru/reference/bm_marginal.md)
  mapper, that is applied before scaling by `weights`.

- A.msk:

  **\[deprecated\]** and has no effect.

- .envir:

  Evaluation environment. Can later be accessed via
  [`bru_comp_env()`](https://inlabru-org.github.io/inlabru/reference/bru_comp_env_extra.md).
  Default: [`parent.frame()`](https://rdrr.io/r/base/sys.parent.html)

- envir_extra:

  Environment for storing
  [`INLA::f()`](https://rdrr.io/pkg/INLA/man/f.html) argument values. If
  NULL, an new environment is created. Can be accessed via
  [`bru_comp_env_extra()`](https://inlabru-org.github.io/inlabru/reference/bru_comp_env_extra.md)

## Value

A `bru_comp` object defining a latent component and its associated
effect mapping.

## Details

As shorthand,
[`bru()`](https://inlabru-org.github.io/inlabru/reference/bru.md) will
understand basic additive formulae describing fixed effect models. For
instance, the components specification `y ~ x` will define the linear
combination of an effect named `x` and an intercept to the response `y`
with respect to the likelihood family stated when calling
[`bru()`](https://inlabru-org.github.io/inlabru/reference/bru.md).
Mathematically, the linear predictor \\\eta\\ would be written as

\$\$\eta = \beta x + c,\$\$

where:

- \\c\\:

  is the *intercept*

- \\x\\:

  is a *covariate*

- \\\beta\\:

  is a *latent variable* associated with \\x\\ and

- \\\psi = \beta x\\:

  is called the *effect* of \\x\\

A problem that arises when using this kind of R formula is that it does
not clearly reflect the mathematical formula. For instance, when
providing the formula to inla, the resulting object will refer to the
random effect \\\psi = \beta \* x \\ as `x`. Hence, it is not clear when
`x` refers to the covariate or the effect of the covariate.

The `bru_comp.character` method is inlabru's equivalent to `INLA`'s
[`f()`](https://rdrr.io/pkg/INLA/man/f.html) function but adds
functionality that is unique to inlabru.

The `main`, `weights`, `group`, `replicate`, and `*_layer` arguments
support `tidy` evaluation expressions, with `.data` and `.env` pronouns.

## Functions

- `bru_component()`: Backwards compatibility alias for `bru_comp()`

## Naming random effects

In INLA, the [`f()`](https://rdrr.io/pkg/INLA/man/f.html) notation is
used to define more complex models, but a simple linear effect model can
also be expressed as

- `formula = y ~ f(x, model = "linear")`,

where [`f()`](https://rdrr.io/pkg/INLA/man/f.html) is the inla specific
function to set up random effects of all kinds. The underlying predictor
would again be \\\eta = \beta \* x + c\\ but the result of fitting the
model would state `x` as the random effect's name. bru allows rewriting
this formula in order to explicitly state the name of the random effect
and the name of the associated covariate. This is achieved by replacing
`f` with an arbitrary name that we wish to assign to the effect, e.g.

- `components = y ~ psi(x, model = "linear")`.

Being able to discriminate between \\x\\ and \\\psi\\ is relevant
because of two functionalities bru offers. The formula parameters of
both [`bru()`](https://inlabru-org.github.io/inlabru/reference/bru.md)
and the prediction method
[predict.bru](https://inlabru-org.github.io/inlabru/reference/predict.md)
are interpreted in the mathematical sense. For instance, `predict` may
be used to analyze the analytical combination of the covariate \\x\\ and
the intercept using

- `predict(fit, data.frame(x=2)), ~ exp(psi + Intercept)`.

which corresponds to the mathematical expression e ^(β + c).

On the other hand, predict may be used to only look at a transformation
of the latent variable \\\beta\_\psi\\

- `predict(fit, NULL, ~ exp(psi_latent))`.

which corresponds to the mathematical expression e ^(β).

## See also

[`bru_input()`](https://inlabru-org.github.io/inlabru/reference/bru_input.md),
[`summary.bru_comp()`](https://inlabru-org.github.io/inlabru/reference/summary.bru_comp.md)

Other component constructors:
[`bru_comp_list()`](https://inlabru-org.github.io/inlabru/reference/bru_comp_list.md)

## Author

Fabian E. Bachl <bachlfab@gmail.com> and Finn Lindgren
<Finn.Lindgren@gmail.com>

## Examples

``` r
# As an example, let us create a linear component. Here, the component is
# called "myLinearEffectOfX" while the covariate the component acts on is
# called "x". Note that a list of components is returned because the
# formula may define multiple components

cmp <- bru_comp_list(~ myLinearEffectOfX(main = x, model = "linear"))
summary(cmp)
#> Label:   myLinearEffectOfX 
#>   Type:  main = linear 
#>   Map:   pipe = multi(main = {autodetect(x)}) 
#>   INLA formula:  
#>     ~ . + f(myLinearEffectOfX, model =
#>       BRU_myLinearEffectOfX_main_model) 
#> Label:   Intercept 
#>   Type:  main = linear 
#>   Map:   pipe = multi(main = {autodetect(1)}) 
#>   INLA formula:  
#>     ~ . + f(Intercept, model = BRU_Intercept_main_model) 
# Equivalent shortcuts:
cmp <- bru_comp_list(~ myLinearEffectOfX(x, model = "linear"))
cmp <- bru_comp_list(~ myLinearEffectOfX(x))
# Individual component
cmp <- bru_comp("myLinearEffectOfX", main = x, model = "linear")
summary(cmp)
#> Label:   myLinearEffectOfX 
#>   Type:  main = linear 
#>   Map:   pipe = multi(main = {autodetect(x)}) 
#>   INLA formula:  
#>     ~ . + f(myLinearEffectOfX, model =
#>       BRU_myLinearEffectOfX_main_model) 
# \donttest{
if (bru_safe_inla()) {
  # As an example, let us create a linear component. Here, the component is
  # called "myEffectOfX" while the covariate the component acts on is called
  # "x":

  cmp <- bru_comp("myEffectOfX", main = x, model = "linear")
  summary(cmp)

  # A more complicated component:
  cmp <- bru_comp("myEffectOfX",
    main = x,
    model = INLA::inla.spde2.matern(fmesher::fm_mesh_1d(1:10))
  )

  # Compound fixed effect component, where x and z are in the input data.
  # The formula will be passed on to MatrixModels::model.Matrix:
  cmp <- bru_comp("eff", ~ -1 + x:z, model = "fixed")
  summary(cmp)
}
#> Label:   eff 
#>   Type:  main = fixed 
#>   Map:   pipe = multi(main = {autodetect(~-1 + x:z)}) 
#>   INLA formula:  
#>     ~ . + f(eff, model = BRU_eff_main_model, hyper =
#>       BRU_eff_main_fixed_hyper) 
# }
```
