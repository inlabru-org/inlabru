# inlabru 2.10.0

## Feature updates

* Add new `ibm_simplify()` generic to handle mapper simplification more generally;
  needed to properly support non-linear component mappers. (version `2.9.0.9004`)
* Add new `bru_mapper_marginal()` mapper class that can be used as part of component
  mapper pipelines. (version `2.9.0.9004`)
* Add new `ibm_eval2()` generic that computes both evaluation and Jacobian,
  avoiding double-computing of the Jacobian, when practical. (version `2.9.0.9005`)
* Add new `bru_timings_plot()` function that plots the time used for each nonlinear iteration
(version `2.9.0.9007`)
* Speed up `bru_fill_missing()` (by orders of magnitude) by changing method for
  finding the nearest available data point. (version `2.9.0.9011`)
* Add new `bru_mapper_shift()` mapper class that works like `bru_mapper_scale()`
  but for additive shifts instead of multiplicative scaling. (version `2.9.0.9012`)
* Added more checks for invalid component or predictor evaluations, to help
  catch user errors sooner, and with more informative messages. (version `2.9.0.9013`)
* Expand `bru_mapper_matrix`, previously used only for component `model = "fixed",
  to allow integer indexing in addition to the previous factor/character-only indexing.
  (version `2.9.0.9014`)

## Bug fixes

* The `is_linear` flag wasn't correctly set for `bru_mapper_logsumexp` mappers.
  Since previous versions did not accept non-linear component mappers, this
  is unlikely to have affected any user code. (Fixed in version `2.9.0.9001`)
* Improved error messages for missing or incomplete LGCP domain specification.
  (version `2.9.0.9002` and `2.9.0.9006`)
* Allow `NULL` in automatic component usage detection. (version `2.9.0.9003`)
* Corrected the crs information for `gorillas$plotsample$counts` and
  `gorillas_sf$plotsample$counts` from `+units=m` to `+units=km`. (version `2.9.0.9010`)
  The geometry information in `counts` is unlikely to have been used in examples
  or analysis code, as the problem would have been immediately obvious;
  plotting or other geometric operations that use the crs information would
  heve been completely wrong, and is only detected now that more code uses the
  crs information at all. Thanks to Dmytro Perepolkin for reporting in issue #205
* Fix problem in `bru_fill_missing()` for cases where the input data object also
  has missing values. (version `2.9.0.9011`)
* Make `eval_spatial()` transform the `where` coordinates to the same crs as the
  input data, for `SpatRaster` and `sf` inputs, to allow different crs specifications.
  (version `2.9.0.9012`)

# inlabru 2.9.0

## Feature updates

* Conversion of code to use `fmesher` for mesh and geometry handling;
  the interface supports existing objects and methods.
  See https://inlabru-org.github.io/fmesher/articles/inla_conversion.html for
  more information.
* General speed improvements, see below for details.
* Added `gg.sf()` method.
* Add experimental support for `stars` via `eval_spatial()`.
  (version `2.8.0.9007`)
* Move the `sp` package from 'Depends' to 'Imports'.  This means that user code
  should either use `sp::` or `library("sp")` to access `sp` methods.
  The `bru_safe_sp()` helper function can be used to check for a safe
  `sp` package configuration during the transition from `rgdal` to `sf`, and
  is only needed if you may run on systems with `sp` installations older than
  "2.0-0" or with `sp::get_evolution_status() < 2`. (version `2.8.2011`)
* Now preserves the previous log output when using `bru_rerun()`,
  and `bru_log()` is now a set of S3 methods, supporting extracting the
  full inlabru log as well `bru`-object specific logs (version `2.8.0.9008`).
  
  Note: From version `2.9.0`, use `bru_log()` to access the global log, and
  `bru_log(fit)` to access a stored estimation log.
  
  Up to version `2.8.0`, `bru_log()` was a deprecated alias for
  `bru_log_message()`. When running on `2.8.0` or earlier, use `bru_log_get()`
  to access the global log, and `cat(fit$bru_iinla$log, sep = "\n")` to print
  a stored estimation object log.

## Bug fixes and speed improvements

* Covariate object component inputs of type `SpatialPolygonsDataFrame`
  were not automatically passed on to `eval_spatial()`. The logic has now changed
  so that any object with a `eval_spatial()` method will trigger a call to
  `eval_spatial()`. See `?input_eval` for further information.
  (version `2.8.0.9001`)
* `fm_crs_is_null()`, `fm_transform()` now supports oblique `fm_crs` CRS objects,
  and `is.na()` methods for the `fm_crs` and `inla.CRS` classes have been added.
  (version `2.8.0.9003`)
* Significant speed up `predict()` by using `quantile(..., names = FALSE)`.
  (version `2.8.0.9004`)
* Improved `row_kron()` code, causing speedups of a factor 2-30 in randomised
  test cases. (version `2.8.0.9005`)
* Removed incorrect code for `sf` method for `eval_spatial()`, causing failure
  when extracting from multiple layers in a single call.
  (version `2.8.0.9007`)
* Improved handling of posterior sample variable extraction in `generate()`
  and `predict()`. Now much faster for large models. (version `2.8.0.9009`)
* Fixed linearisation issue when using only the `*_latent` form of a component.
  (version `2.8.0.9015`)
* Workaround for equivalent but textually different CRS/WKT information in
  `bru_fill_missing()`. (version `2.8.0.9016`, fixes #200)

## Deprecation of old functions

* `eval_SpatialDF` removed, deprecated since `2.8.0`. See `eval_spatial` instead.
* `stransform`, `ibm_amatrix`, `ibm_valid_input` removed, deprecated since `2.7.0`.
  See `fm_transform` and `ibm_jacobian` instead.
* `bru_mapper_offset`, deprecated since `2.6.0` now returns a pure `bru_mapper_const`
  object, and all `bru_mapper_offset` `ibm_*` methods have been removed.
* `init.tutorial` removed, deprecated since `2.5.0`
* `generate.inla` and `predict.inla` removed, deprecated since `2.1.0`

# inlabru 2.8.0

## Feature updates

* The iterative inla method has been given both sharper internal `inla()` optimisation
  criteria for the iterations (thanks to Haavard Rue), _and_ a more relaxed
  nonlinear iteration stopping criterion; the default `bru_method$rel_tol`
  values has been changed from 1 to 10 percent change. The iterations are
  terminated when all latent and hyper-parameter mode changes fullfil
  `|change|/SD < rel_tol`, and the non-linear line search is inactive.
  This seems to strike a useful balance between the different optimisation
  criteria, allowing the iterations to converge faster and also detect that
  convergence sooner.
* The logic for which components are needed for a predictor expression
  (in `like()` or `generate()`/`predict()`) has been updated to when possible
  extract the list of components from the expression itself.
  The user can override this default if necessary, using the `include`/`exclude` arguments.
  
  The `bru_used()` methods are used to guess the needed component names, applied
  to the right-hand side of the `formula` arguments.  The `allow_latent` argument
  to `like()` has been deprecated in favour of `include_latent`
  (by default auto-detected for use of `_latent` and `_eval`).
  
  The internal information storage is handled by the new `bru_used()`
  methods, that can also be used directly by the user and supplied via the
  `used` argument to `like()`/`generate()`/`predict()`.
* Add `fm_int()` integration methods, replacing the old `ipmaker()` and `ipoints()` methods.
  Supports both `sf` and `sp` sampler objects.
* Add `fm_pixels()` methods for gridded points. The old 
  `pixels()` method now calls `fm_pixels(..., format = "sp")`
* `eval_spatial` support for sf objects (for point-in-polygon data lookups)
* Allow precomputed spatial covariates in the data for point process observations
* Add `edge|int|ext.linewidth` arguments to `gg.inla.mesh` #188
* Rename the `predict()` and `generate()` `data` arguments to `newdata`, for
  better compatibility with other `predict()` methods.  The old argument name
  will still be accepted, but give a warning.  Code that does not name the `data`
  argument is not affected.
* Note: Coordinate names for `Spatial*` objects have been inconsistently
  available in the predictor expression evaluation. However, due to how internal
  conversions might inadvertently change these names, they can not be relied
  on, and they are no longer being made available to the predictor expression.
  As a side effect, this change also speeds up some `bru()` runs by around a
  factor 2, since it avoids converting the `Spatial*` to a regular `data.frame`
  in time-sensitive core evaluation code.
  
  If you need access to the raw coordinate values, use explicit calls to
  `sp::coordinates(.data.)` (e.g. for custom spatial covariate evaluation.).
  When possible, use the built-in covariate evaluation method, `eval_spatial()`,
  either implicitly with `comp(covariate, ...)` or explicitly,
  `comp(eval_spatial(covariate, where = .data.), ...)`, that handles `crs` information
  correctly.  Also consider transitioning from `sp` to `sf` data storage, using
  `geometry` instead of raw coordinates.

## Bug and dependency updates

* Remove `rgdal` and `maptools` dependencies #178
* Add `bru_safe_sp()` to check if `sp` can be used safely (checks `rgdal` availability
  and `sp` evolution status, optionally forcing use of `sf`) #178
* Remove PROJ4 support #178
* Change `rgl.*` functions to `*3d`. Thanks to Duncan Murdoch #181
* Speed up `ibm_jacobian.bru_mapper_harmonics` for large models
* Add workarounds for inconsistent polygon orientation resulting from `sf::st_*`
  calls that don't account for the `geos` canonical representation being CW,
  whereas the canonical Simple Features representation being CCW. See
  https://github.com/r-spatial/sf/issues/2096
  
# inlabru 2.7.0

## Feature overview

* Added support for `sf ` and `terra` inputs to most methods
* Expanded geometry and mesh handling methods
* Expanded `bru_mapper()` system
* Added convergence diagnostics plot with `bru_convergence_plot()`

## Feature details

* Allow `NA` input for default 1D mappers to generate effect zero, like
  in `inla()`.
* New and expanded methods `fm_crs()`, `fm_CRS()`, `fm_transform()`,
  `fm_ellipsoid_radius()`, and `fm_length_unit()` to further support `sf` objects.
  The `fm_crs()` extraction method also supports `terra` objects.
* `bru_fill_missing()` now supports `terra` `SpatRaster` data and and `sf` locations.
* New experimental methods `fm_evaluator()` and `fm_evaluate()`, replacing the
  `INLA` `inla.mesh.projector` and `inla.mesh.project` methods.
* Experimental integration support for sphere and globe meshes.
* Allow `sf` input to `family="cp"` models.
* Further `bru_mapper()` method updates;

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
  * `summary` methods for `bru_mapper` objects (`summary.bru_mapper()`)
  * Removed `methods` argument from `bru_mapper_define()`.  Implementations
    should register S3 methods instead.
    
## Bug fixes

* Remove unused `spatstat.core` dependency. Fixes #165
* Fixed issue with plain mapper evaluation in the `ibm_eval.default()`
  and `ibm_eval.bru_mapper_collect()` methods, where they would return zeros
  instead of the intended values.
  The main component evaluation and estimation code was not directly affected
  as that is based on the `bru_mapper_multi()` class methods that rely on the
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

