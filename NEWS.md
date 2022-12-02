# inlabru 2.7.0

## Feature overview

* Added support for `sf ` and `terra` inputs to most methods

* Expanded geometry and mesh handling methods

* Expanded `bru_mapper` system

* Added convergence diagnostics plot with `bru_convergence_plot()`

## Feature details

* Allow `NA` input for default 1D mappers to generate effect zero, like
  in `inla()`.
  
* New and expanded methods `fm_crs`, `fm_CRS`, `fm_transform`,
  `fm_ellipsoid_radius`, and `fm_length_unit` to further support `sf` objects.
  The `fm_crs` extraction method also supports `terra` objects.

* `bru_fill_missing` now supports `terra` `SpatRaster` data and and `sf` locations.
  
* New experimental methods `fm_evaluator` and `fm_evaluate`, replacing the
  `INLA` `inla.mesh.projector` and `inla.mesh.project` methods.
  
* Experimental integration support for sphere and globe meshes.
  
* Allow `sf` input to `family="cp"` models.

* Further `bru_mapper` method updates;

  * Deprecated `ibm_amatrix()` and `names()`
  methods, replaced by `ibm_jacobian()` and `ibm_names()`.
  
  * Introduced `bru_mapper_pipe()`, used to link mappers in sequence.
  
  * Introduced `bru_mapper_aggregate()` and `bru_mapper_logsumexp()`,
    used for blockwise weighted sums and log-sum-exp mappings,
    `output[k] = sum(weights[block==k]*state[block==k])))` and
    `output[k] = log(sum(weights[block==k]*exp(state[block==k])))`,
    with optional weight normalisation within each block.  Allows
    providing the weights as log-weights, and uses block-wise shifts to
    avoid potential overflow.
  
  * `summary` methods for `bru_mapper` objects
  
  * Removed `methods` argument from `bru_mapper_define()`.  Implementations
    should register S3 methods instead.
    
## Bug fixes

* Remove unused `spatstat.core` dependency. Fixes #165

* Fixed issue with plain mapper evaluation in the `ibm_eval.default`
  and `ibm_eval.bru_mapper_collect` methods, where they would return zeros
  instead of the intended values.
  The main component evaluation and estimation code was not directly affected
  as that is based on the `bru_mapper_multi` class methods that rely on the
  Jacobians instead.  The bug would therefore mainly have impacted the future,
  not yet supported nonlinear mapper extensions.
  
* Fix for `eval_spatial.SpatRaster`; Work around inconsistent logic in
  `terra::extract(..., layer)` when `length(layer)==1` or `nrow(where)==1`.
  Fixes #169

  * Add `indexed` logical option to `bru_mapper_factor()`, to allow
  factor inputs to be mapped to index values, as needed for `group` and
  `replicate`. Fixes #174

# inlabru 2.6.0

## Features
  
* Add `bru_get_mapper` generic, and associated methods for `inla.spde` and
  `inla.rgeneric` objects. This allows `inlabru` to automatically extract
  the appropriate `bru_mapper` object for each model component, and can be used
  as a hook by external packages implementing new INLA object classes.
  
* Add a `weights` argument for `like()`, for likelihood-specific log-likelihood
  weights, passed on to the `INLA::inla()` weights argument. Evaluated in the
  data context.
  
* The `<component>_eval()` methods available in predictor expressions
  now handle optional scaling weights, like in ordinary component effect
  evaluation.

* Add `terra` support for covariate inputs

* The component `*_layer` arguments are now evaluated in the data context,
  to allow dynamic layer selection for spatial raster covariates.  A new
  generic `eval_spatial()` provides support for grid/pixel based
  `Spatial*DataFrame` evaluation, and `SpatRaster`. Expanded support
  is in progress.

* New vignettes on the `bru_mapper` system, `component` definitions,
  and `prediction_scores`
  
* General overhaul of the `bru_mapper` and linearised predictor system,
  to prepare for new features.
  
  * Add `ibm_eval` generic for evaluating mappers for given states.
  
  * Add `bru_mapper_taylor`, used as an internal mapper for linearised
  mappers. This and `ibm_eval` is aimed at future support for nonlinear
  mappers. Associated new generic methods: `ibm_{is_linear,jacobian,linear}`.
  
  * New mapper implementations should use `ibm_jacobian` instead of `ibm_amatrix`.
    This allows defining a linearised mapper via
   `ibm_eval(input, state0) + ibm_jacobian(input, state0) %*% (state - state0)`.
  
  * New mapper class `bru_mapper_const`, which replaces `bru_mapper_offset`.
  `bru_mapper_offset` is now deprecated and will produce warnings.
  
## Bug fixes

* Enable both datum/ensemble container for ellipsoid information, to support
  `epsg:4326`. Fixes #154
  
* Make duplicated component names an error instead of a warning.
  Relates to #155

* Fix `Tsparse` assumptions in `row_kron` to prepare for Matrix `1.5-2`.
  Fixes #162
  
# inlabru 2.5.3

## Features

* Add `bru_mapper_harmonics` mapper for `cos` and `sin` basis sets.

* Allow `predict()` input data to be be a list.

* Allow arbitrary quantile summaries in `predict()`

* Remove `cv`, `var`, `smin`, `smax` summaries from `predict()`

* Add `mean.mc_std_err` and `sd.mc_std_err` output to `predict()`

* Add `robins_subset` data set and associated variable coefficient web vignette

## Bug fixes

* Propagate multi-likelihood A-matrix information instead of recomputing.
  Fixes iteration issue for bym2 and other `bru_mapper_collect` models.
  
* Turn on predictor summaries during iterations to allow `inla.mode="classic"`
  to use proper line search.
  
* Avoid deprecated Matrix (>=1.4-2) class coercion methods

* Work around for lack of full Matrix and ModelMatrix support for the `unique`
  method. Fixes #145

# inlabru 2.5.2

* More robust package checks

* More robust namespace and INLA availability checks

* Add package vignette with links to the website examples

# inlabru 2.5.1

* Revert to R language features compatible with R 4.0.5

* Use `strategy="gaussian"` during iterations.

# inlabru 2.5.0

## Features

* Add `bru()` timing information in `$bru_timings` and `$bru_iinla$timings`

* Add `SpatialPolygonsDataFrame` support to `gg()` methods

* Allow accessing `E` and `Ntrials` from `response_data` and `data`
  (further special arguments remain to be added)
  
* `deltaIC` improvements

* New transformation helper tools `bru_{forward/inverse}_transformation()`

* Experimental support for matrix and formula component inputs. E.g. with
    `~ name(~ -1 + a + b + a:b, model = "fixed")`, covariate fixed effect interaction
  specifications can be made. For formula input, `MatrixModels::model.Matrix()`
  is called to construct matrix input that is then used as the A-matrix for
  fixed effects, one per column, added up to form the combined effect.

* Documentation and examples improvements

## Bug fixes

* Fix A-matrix construction for `evaluate_model()` for cases where the `inla_f`
  argument matters

* More efficient and robust mesh integration code

* Cleanup of environment handling for component lists

# inlabru 2.4.0

## Features

* Allow predictors to have different size than the input data.  The `data` argument
  is now allowed to be a `list()`, and the new argument `response_data` allows separate
  specification of component inputs and response variables.
  
* Add `bru_mapper_collect` class for handling sequential collections of
  mappers, including collections where all but the first mapper is hidden from the
  `INLA::f()` arguments `n` and `values`, as needed to support e.g. "bym2" models.
  
* Add `control.family` as a direct argument to `like()`. Gives a warning if a
  `control.family` argument is supplied to the the `options` argument of `bru()`,
  but at least one likelihood has `control.family` information. (Issue #109)
  
## Bugfixes

* Fix support for `SpatialPointsDataFrame` and `SpatialGridDataFrame` input
  to `bru_fill_missing()`

* Force explicit `model = "offset"` components instead of special options, to
  avoid interfering with the linearisation system (Issue #123)
  
* Make the iterations more robust by resetting the internal INLA predictor
  states to initial value zero at each step
  
## Miscellaneous
  
* Rename the option `bru_method$stop_at_max_rel_deviation` to `bru_method$rel_tol`.
  Automatic conversion to the new name, but a warning is given.

* Add option `bru_method$max_step` to control the largest allowed line search
  scaling factor. See `?bru_options`
  
* New default option `bru_compress_cp` set to `TRUE` to compress the predictor
  expression for `family="cp"` to use a single element for the linear predictor sum.
  
# inlabru 2.3.1

* Documentation and dependency updates for CRAN compatibility

* See NEWS for version 2.3.0 for the major updates since version 2.1.13

# inlabru 2.3.0

## Breaking changes since version 2.1.13

* The model component argument `map` has been deprecated. Use `main` to specify
  the main component input, `~ elev(main = elevation, model = "rw2")`.
  Unlike the old `map` argument, `main` is the first one, so the shorter version
  `~ elev(elevation, model = "rw2")` also works.

* Intercept-like components should now have explicit inputs, e.g. `~ Intercept(1)`
  to avoid accidental confusion with other variables.
  
* The argument list for `bru()` has been simplified, so that all arguments except
  `components` and `options` must either be outputs from calls to `like()`, or
  arguments that can be sent to a single `like()` call.
  
* The option setting system has been replaced with a more coherent system;
  see `?bru_options()` for details.
  
* The `samplers` and `domain` system for `lgcp` models is now stricter, and
  requires explicit `domain` definitions for all the point process dimensions.
  Alternatively, user-defined integration schemes can be supplied via the `ips`
  argument.

## New features since version 2.1.13

* The model component input arguments `main`, `group`, `replicate`, and `weights`
  can now take general R expressions using the data inputs. Special cases are detected:
  `SpatialPixels/GridDataFrame` objects are evaluated at spatial locations if
  the input data is a `SpatialPointsDataFrame` object. Functions are evaluated
  on the data object, e.g. `field(coordinates, model = spde)`

* The component arguments `mapper`, `group_mapper`, and `replicate_mapper` can be
  used for precise control of the mapping between inputs and latent variables.
  See `?bru_mapper` for more details. Mapper information is automatically extracted
  from `INLA::inla.spde2.pcmatern()` model objects.
  
* The R-INLA `weights` and `copy` features are now supported.

* The predictor expressions can access the data object directly via `.data.`

* If data from several rows can affect the same output row, the `allow_combine = TRUE`
  argument must be supplied to `like()`

* The `include` and `exclude` arguments to `like()`, `generate()`, and `predict()`
  can be used to specify which components are used for a given likelihood model
  or predictor expression. This can be used to prevent evaluation of components
  that are invalid for a likelihood or predictor.
  
* Predictor expressions can access the latent state of a model component directly,
  by adding the suffix `_latent` to the component name, e.g. `name_latent`.
  For `like()`, this requires
  `allow_latent = TRUE` to activate the needed linearisation code for this.
  
* Predictor expressions can evaluate component effects for arbitrary inputs by
  adding the suffix `_eval` to access special evaluator functions, e.g.
  `name_eval(1:10)`. This is useful for evaluating the 1D effect of spatial covariates.
  See the NEWS item for version 2.2.8 for further details.

* The internal system for predictor linearisation and iterated INLA inference
  has been rewritten to be faster and more robust
  
* See the NEWS entries for versions 2.1.14 to 2.2.8 for further details on new
  features and bug fixes

# inlabru 2.2.8

* Add `_eval` suffix feature for `generate.bru` and `predict.bru`, that
  provides a general evaluator function for each component, allowing evaluation
  of e.g. nonlinear effects of spatial covariates as a function of the covariate
  value instead of the by the spatial evaluator used in the component definition.
  For example, with `components = ~ covar(spatial_grid_df, model = "rw1")`, the
  prediction expression can have `~ covar_eval(covariate)`, where `covariate`
  is a data column in the prediction data object.
  
  For components with `group` and `replicate` features, these also need to be
  provided to the `_eval` function, with
  `..._eval(..., group = ..., replicate = ...)`
  
  This feature is built on top of the `_latent` suffix feature, that gives
  direct access to the latent state variables of a component, so in order to
  use `_eval` in the model predictor itself, you must use
  `like(..., allow_latent = TRUE)` in the model definition.

# inlabru 2.2.7

* Add support for `ngroup` and `nrep` in component definitions

* Updated `mexdolphin` and `mrsea` data sets, with consistent km units and
  improved mesh designs

# inlabru 2.2.6

* Add `predict(..., include)` discussion to distance sampling vignette, for
  handling non-spatial prediction in spatial models.
  
* Fix bugs in `gg.SpatialLines`

# inlabru 2.2.5

* Vignette corrections

* Documentation improvements

* Fix minor bug in `Spatial*` object handling and plotting

# inlabru 2.2.4

* Properly extract the joint latent conditional mode instead of the
  marginal latent conditional mode

# inlabru 2.2.2

* Fixed issue with `predict()` logic for converting output to `Spatial*DataFrame`

* Use `control.mode=list(restart=FALSE)` in the final inla run for nonlinear
  models, to avoid an unnecessary optimisation.

* Fix issues in `pixels()` and `bru_fill_missing()` for `Spatial*DataFrame`
  objects with `ncol=0` data frame parts.

# inlabru 2.2.1

* Fixed code regression bug for function input of covariates

# inlabru 2.2.0

* Support for the INLA "copy" feature, `comp2(input, copy = "comp1")`

* Allow component weights to be an unnamed parameter, `comp(input, weights, ...)`

* Direct access to the data objects in component inputs and predictor
  expressions, as `.data.`, allowing e.g. `covar(fun(.data.), ...)` for a complex
  covariate extractor method `fun()`
  
* Partial support for spherical manifold meshes

* Uses INLA integration strategy "eb" for initial nonlinear iterations, and a
  specified integration strategy only for the final iteration, so that the
  computations are faster, and uses the conditional latent mode as
  linearisation point.

# inlabru 2.1.15

* New options system

* New faster linearisation method

* New line search method to make the nonlinear inla iterations robust

* Method for updating old stored estimation objects

* System for supplying mappings between latent models and evaluated effects
  via `bru_mapper` objects

* Improved factor support; Either as "contrast with the 1st level", via the
  special `"factor_contrast"` model, or all levels with model `"factor_full"`.
  Further options planned (e.g. a simpler options to fix the precision
  parameter).  The estimated coefficients appear as random effects in the
  `inla()` output.

* Interface restructuring to support new features while keeping most
  backwards compatibility. Change `map=` to `main=` or unnamed first argument;
  Since `main` is the first parameter, it doesn't need to be a named argument.

* Keep components with zero derivative in the linearisation

* PROJ6 support

* Add random seed option for posterior sampling

* Add package unit testing

* New backend code to make extended feature support easier

* New `int.args` option to control spatial integration resolution,
  thanks to Martin Jullum (`martinju`)

# inlabru 2.1.13

* Fix CRAN complaint regarding documentation

# inlabru 2.1.12

* Workaround an integration points error for old (ca pre-2018) INLA versions

# inlabru 2.1.11

* Add CITATION file

# inlabru 2.1.10

* Fix internal sampling bug due to INLA changes

# inlabru 2.1.9

* Remove unused `VignetteBuilder` entry from `DESCRIPTION`

# inlabru 2.1.8

* Update default options

* Prevent `int.polygon` from integrating outside the mesh domain,
  and generally more robust integration scheme construction.

* Fix `bru()` to `like()` parameter logic. (Thanks to Peter Vesk for bug example)

# inlabru 2.1.7

* Added a `NEWS.md` file to track changes to the package.

* Added `inla` methods for `predict()` and `generate()` that convert
  `inla` output into `bru` objects before calling the `bru` prediction
  and posterior sample generator.

* Added protection for examples requiring optional packages

* Fix `sample.lgcp` output formatting, extended CRS support, and more efficient sampling algorithm

* Avoid dense matrices for effect mapping

# inlabru 2.1.4

* `iinla()` tracks convergence of both fixed and random effects

# inlabru 2.1.3

* Added matrix geom `gg.matrix()`

* Fixed CRAN test issues

