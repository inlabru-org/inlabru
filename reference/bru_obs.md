# Observation model construction for [`bru()`](https://inlabru-org.github.io/inlabru/reference/bru.md)

Observation model construction for usage with
[`bru()`](https://inlabru-org.github.io/inlabru/reference/bru.md).

Note: Prior to version `2.12.0`, this function was called `like()`, and
that alias will remain for a while until examples etc have been updated
and users made aware of the change. The name change is to avoid issues
with namespace clashes, e.g. with
[`data.table::like()`](https://rdrr.io/pkg/data.table/man/like.html),
and also to signal that the function defines observation models, not
just likelihood functions.

## Usage

``` r
bru_obs(
  formula = . ~ .,
  family = "gaussian",
  data = NULL,
  response_data = NULL,
  data_extra = NULL,
  E = NULL,
  Ntrials = NULL,
  weights = NULL,
  scale = NULL,
  domain = NULL,
  samplers = NULL,
  ips = NULL,
  used = NULL,
  is_rowwise = NULL,
  aggregate = NULL,
  aggregate_input = NULL,
  response_block = NULL,
  control.family = NULL,
  control.gcpo = NULL,
  tag = NULL,
  options = list(),
  .envir = parent.frame(),
  include = deprecated(),
  exclude = deprecated(),
  include_latent = deprecated(),
  allow_combine = deprecated()
)

like(
  ...,
  E = NULL,
  Ntrials = NULL,
  weights = NULL,
  scale = NULL,
  options = list(),
  .envir = parent.frame()
)

bru_obs_list(..., .tag = NULL)

# S3 method for class 'list'
bru_obs_list(object, ..., .tag = NULL)

# S3 method for class 'bru_obs'
bru_obs_list(..., .tag = NULL)

# S3 method for class 'bru_obs_list'
bru_obs_list(..., .tag = NULL)

# S3 method for class 'bru_obs'
c(...)

# S3 method for class 'bru_obs_list'
c(...)

# S3 method for class 'bru_obs_list'
x[i]

like_list(...)

bru_like_list(...)
```

## Arguments

- formula:

  a `formula` where the right hand side is a general R expression
  defines the predictor used in the model.

- family:

  A string identifying a valid
  [`INLA::inla`](https://rdrr.io/pkg/INLA/man/inla.html) likelihood
  family. The default is `gaussian` with identity link. Apart from the
  likelihoods provided by INLA (see
  `names(INLA::inla.models()$likelihood)`) inlabru supports

  `cp`

  :   Cox process likelihood, for fitting point process models. The
      [`lgcp()`](https://inlabru-org.github.io/inlabru/reference/lgcp.md)
      function is a shortcut to `bru(..., family = "cp")`.

  `nzbinomial`,`nzbetabinomial`,`nznbinomial`,`nzcenpoisson`

  :   Non-zero versions of "binomial", "betabinomial", "nbinomial", and
      "cenpoisson", respectively, implemented by setting the zero
      probability parameter close to zero in the corresponding
      "zeroinflated\*0" models. It will check if these models have
      native INLA implementations, and use those instead if available.
      In INLA `26.06.07`, only "nzpoisson" has a native implementation.

- data:

  Predictor expression-specific data, as a `data.frame`, `tibble`, or
  `sf`. Since `2.12.0.9023`, deprecated support for
  `SpatialPoints[DataFrame]` objects.

- response_data:

  Observation/response-specific data for models that need different
  size/format for inputs and response variables, as a `data.frame`,
  `tibble`, or `sf`. Since `2.12.0.9023`, deprecated support for
  `SpatialPoints[DataFrame]` objects.

- data_extra:

  object convertible with
  [`as.list()`](https://rdrr.io/r/base/list.html) with additional
  variables to be made available in predictor evaluations. Variables
  with the same names as the data object will be ignored, unless
  accessed via `.data_extra.[["name"]]` or `.data_extra.$name` in the
  formula.

- E:

  Exposure/effort parameter for family = 'poisson' passed on to
  [`INLA::inla`](https://rdrr.io/pkg/INLA/man/inla.html). Special case
  if family is 'cp': rescale all integration weights by a scalar `E`.
  For sampler specific reweighting/effort, use a `weight` column in the
  `samplers` object instead, see
  [`fmesher::fm_int()`](https://inlabru-org.github.io/fmesher/reference/fm_int.html).
  Default taken from `options$E`, normally `1`.

- Ntrials:

  A vector containing the number of trials for the 'binomial'
  likelihood. Default taken from `options$Ntrials`, normally `1`.

- weights:

  Fixed (optional) weights parameters of the likelihood, so the
  log-likelihood`[i]` is changed into `weights[i] * log_likelihood[i]`.
  Default value is `1`. WARNING: The normalizing constant for the
  likelihood is NOT recomputed, so ALL marginals (and the marginal
  likelihood) must be interpreted with great care.

  For `family = "cp"`, the weights are applied as `sum(weights * eta)`
  in the point location contribution part of the log-likelihood, where
  `eta` is the linear predictor, and do not affect the integration part
  of the likelihood. This can be used to implement approximative methods
  for point location uncertainty.

- scale:

  Fixed (optional) scale parameters of the precision for several models,
  such as Gaussian and student-t response models.

- domain, samplers, ips:

  Arguments used for `family="cp"` and `aggregate=`.

  `domain`

  :   Named list of domain definitions, see
      [`fmesher::fm_int()`](https://inlabru-org.github.io/fmesher/reference/fm_int.html).

  `samplers`

  :   Integration domain for `family="cp"` or subdomains for
      `aggregate=`, see
      [`fmesher::fm_int()`](https://inlabru-org.github.io/fmesher/reference/fm_int.html).

  `ips`

  :   Integration points. Defaults to
      [fmesher::fm_int](https://inlabru-org.github.io/fmesher/reference/fm_int.html)`(domain, samplers)`.
      If explicitly given, overrides `domain` and `samplers`.
      [`fmesher::new_fm_int()`](https://inlabru-org.github.io/fmesher/reference/new_fm_int.html)
      (from fmesher `0.5.0.9013`) can be used for manually constructed
      integration schemes.

- used:

  Either `NULL` (default) or a
  [`bru_used()`](https://inlabru-org.github.io/inlabru/reference/bru_used.md)
  object. When, `NULL`, the information about what effects and latent
  vectors are made available to the predictor evaluation is defined by
  `bru_used(formula)`, which will include all effects and latent vectors
  used by the predictor expression.

- is_rowwise, allow_combine:

  logical; If `is_rowwise` is `FALSE`, the predictor expression may
  involve several rows of the input data to influence the same row. When
  `NULL`, it defaults to `TRUE`, unless `response_data` is non-`NULL`,
  or `data` is a `list`, or the likelihood construction requires it. The
  `allow_combine` argument is a deprecated and inverted version, such
  that `is_rowwise = !allow_combine`.

- aggregate:

  character ("none" or a valid name for the `type` argument of
  [`bm_aggregate()`](https://inlabru-org.github.io/inlabru/reference/bm_aggregate.md))
  or an aggregation `bru_mapper` object
  ([`bm_aggregate()`](https://inlabru-org.github.io/inlabru/reference/bm_aggregate.md),
  [`bm_logsumexp()`](https://inlabru-org.github.io/inlabru/reference/bm_logsumexp.md),
  or
  [`bm_logitaverage()`](https://inlabru-org.github.io/inlabru/reference/bm_logitaverage.md)).
  Default `NULL`, interpreted as "none". **\[experimental\]**, available
  from version `2.12.0.9013`.

- aggregate_input, response_block:

  `NULL` or an optional input list to the mapper defined by non-NULL
  `aggregate`, overriding the default,

      aggregate_input = list(
        block = .data.[[".block"]],
        weights = .data.[["weight"]],
        n_block = bru_response_size(.response_data.)
      ),
      response_block = .response_data.[[".block"]]

  **\[experimental\]**, available from version `2.12.0.9013`.

  From `2.13.0.9016` to `2.14.1.9007`, it would look for a
  `block_response` character element in the list, but from
  `2.14.1.9008`, the separate argument `reponse_block` should be used
  instead, with an expression to be evaluated in the input data context.
  `response_block` should evaluate to a vector of the same length as the
  response data, with values that can be used to index into the rows of
  the response variable to determine which rows of the predictor
  expression to aggregate, allowing character or factor aggregation
  block information to be used by replacing `block` with
  `match(block, response_block)`. If not supplied, the variable `.block`
  in the `response_data` is tried. If that isn't available, the `block`
  information must be supplied directly as a numeric or integer vector,
  indexing into the rows of the response variable. Having no
  `response_block` variable is equivalent to having
  `response_data$.block = seq_len(bru_response_size(response_data))`.

- control.family:

  A optional `list` of
  [`INLA::control.family`](https://rdrr.io/pkg/INLA/man/control.family.html)
  options

- control.gcpo:

  A optional `list` of
  [`INLA::control.gcpo`](https://rdrr.io/pkg/INLA/man/control.gcpo.html)
  options

- tag:

  character; Name that can be used to identify the relevant parts of
  INLA predictor vector output, via
  [`bru_index()`](https://inlabru-org.github.io/inlabru/reference/bru_index.md).

- options:

  A
  [bru_options](https://inlabru-org.github.io/inlabru/reference/bru_options.md)
  options object or a list of options passed on to
  [`bru_options()`](https://inlabru-org.github.io/inlabru/reference/bru_options.md)

- .envir:

  The evaluation environment to use for special arguments (`E`,
  `Ntrials`, `weights`, and `scale`) if not found in `response_data` or
  `data`. Defaults to the calling environment.

- include, exclude, include_latent:

  **\[deprecated\]**, use `used` instead.

- ...:

  For `bru_obs_list.bru_obs`, one or more `bru_obs` objects

- .tag:

  Optional name to assign to a single `bru_obs` object. Reserved for
  internal use.

- object:

  A list of `bru_obs` and/or `bru_obs_list` objects

- x:

  `bru_obs_list` object from which to extract element(s)

- i:

  indices specifying elements to extract

## Value

A likelihood configuration which can be used to parameterise
[`bru()`](https://inlabru-org.github.io/inlabru/reference/bru.md).

## Details

The `E`, `Ntrials`, `weights`, and `scale` arguments are evaluated in
the data context, with values from `response_data` taking precedence
over `data`.

## Methods (by generic)

- `bru_obs_list(bru_obs)`: Combine one or more lists of `bru_obs`
  observation model objects into a `bru_obs_list` object

- `c(bru_obs)`: Combine several `bru_obs` objects into a `bru_obs_list`
  object

## Functions

- `like()`: **\[deprecated\]** Legacy `like()` method for `inlabru`
  prior to version `2.12.0`. Use `bru_obs()` instead.

- `bru_obs_list()`: Combine `bru_obs` observation model object into a
  `bru_obs_list` object

- `bru_obs_list(list)`: Combine one or more lists of `bru_obs`
  observation model objects into a `bru_obs_list` object

- `bru_obs_list(bru_obs_list)`: Combine one or more `bru_obs_list`
  objects into a `bru_obs_list` object

- `c(bru_obs_list)`: Combine several `bru_obs_list` objects into a
  `bru_obs_list` object

- `like_list()`: **\[deprecated\]** Backwards compatibility for versions
  `<= 2.12.0`. For later versions, use
  [`as_bru_obs_list()`](https://inlabru-org.github.io/inlabru/reference/as_bru_obs.md),
  `bru_obs_list()`, or [`c()`](https://rdrr.io/r/base/c.html).

- `bru_like_list()`: **\[deprecated\]** Backwards compatibility for
  versions `<= 2.12.0.9017`. For later versions, use
  [`as_bru_obs_list()`](https://inlabru-org.github.io/inlabru/reference/as_bru_obs.md),
  `bru_obs_list()` or [`c()`](https://rdrr.io/r/base/c.html).

## See also

[`bru_response_size()`](https://inlabru-org.github.io/inlabru/reference/bru_response_size.md),
[`bru_used()`](https://inlabru-org.github.io/inlabru/reference/bru_used.md),
[`bru_comp()`](https://inlabru-org.github.io/inlabru/reference/bru_comp.md),
[`bru_comp_eval()`](https://inlabru-org.github.io/inlabru/reference/bru_comp_eval.md)

[`summary.bru_obs()`](https://inlabru-org.github.io/inlabru/reference/bru_obs_print.md)

## Author

Fabian E. Bachl <bachlfab@gmail.com>

Finn Lindgren <finn.lindgren@gmail.com>

## Examples

``` r
# \donttest{
if (bru_safe_inla() &&
  require(ggplot2, quietly = TRUE) &&
  require(patchwork, quietly = TRUE)) {
  # The 'bru_obs()' (previously 'like()') function's main purpose is to set
  # up observation models, both for single- and multi-likelihood models.
  # The following example generates some random covariates which are observed
  # through two different random effect models with different likelihoods

  # Generate the data

  set.seed(123)

  n1 <- 200
  n2 <- 10

  x1 <- runif(n1)
  x2 <- runif(n2)
  z2 <- runif(n2)

  y1 <- rnorm(n1, mean = 2 * x1 + 3)
  y2 <- rpois(n2, lambda = exp(2 * x2 + z2 + 3))

  df1 <- data.frame(y = y1, x = x1)
  df2 <- data.frame(y = y2, x = x2, z = z2)

  # Single likelihood models and inference using bru are done via

  cmp1 <- y ~ -1 + Intercept(1) + x
  fit1 <- bru(cmp1, family = "gaussian", data = df1)
  summary(fit1)

  cmp2 <- y ~ -1 + Intercept(1) + x + z
  fit2 <- bru(cmp2, family = "poisson", data = df2)
  summary(fit2)

  # A joint model has two likelihoods, which are set up using the bru_obs
  # function

  lik1 <- bru_obs(
    "gaussian",
    formula = y ~ x + Intercept,
    data = df1,
    tag = "norm"
  )
  lik2 <- bru_obs(
    "poisson",
    formula = y ~ x + z + Intercept,
    data = df2,
    tag = "pois"
  )

  # The union of effects of both models gives the components needed to run
  # bru

  jcmp <- ~ x + z + Intercept(1)
  jfit <- bru(jcmp, lik1, lik2)

  bru_index(jfit, "norm")
  bru_index(jfit, "pois")

  # Compare the estimates

  p1 <- ggplot() +
    gg(fit1$summary.fixed, bar = TRUE) +
    ylim(0, 4) +
    ggtitle("Model 1")
  p2 <- ggplot() +
    gg(fit2$summary.fixed, bar = TRUE) +
    ylim(0, 4) +
    ggtitle("Model 2")
  pj <- ggplot() +
    gg(jfit$summary.fixed, bar = TRUE) +
    ylim(0, 4) +
    ggtitle("Joint model")

  (p1 / p2 / pj)

  # Non-zero binomial example:
  nzdata <- data.frame(ntrials = rep(2:6, 10))
  nzdata$x <- rnorm(nrow(nzdata))
  nzdata$count <- rbinom(
    nrow(nzdata),
    size = nzdata$ntrials,
    prob = plogis(nzdata$x)
  )
  nzdata <- nzdata[nzdata$count > 0, ]
  truncated_binomial_obs <-
    bru_obs(
      formula = count ~ x,
      family = "nzbinomial",
      data = nzdata,
      Ntrials = ntrials
    )
  fit <- bru(~ 0 + x, truncated_binomial_obs)
}
# }
```
