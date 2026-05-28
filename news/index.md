# Changelog

## inlabru 2.14.1

CRAN release: 2026-05-01

### Bug fixes

- Fix issue [\#293](https://github.com/inlabru-org/inlabru/issues/293),
  where automated mapper construction for `factor_contrast` always used
  the alphabetically first level as contrast, instead of the
  factor-defined first level (`2.14.0.9001`)
- Fix `eval_spatial<SpatRaster>` to return factor data for pure factor
  layer extraction. Layer-specific
  [`terra::extract()`](https://rspatial.github.io/terra/reference/extract.html)
  may return integer instead of factor data, as multi-layer extraction
  would otherwise need to mix data of different types (`2.14.0.9002`)
- Check for “error” inheritance for input evaluation and other errors
  instead of “simpleError”, as R \> 4.6 uses more fine-grained error
  sub-classes in more cases than before (`2.14.0.9004`)

### General updates

- Move `fmesher` from `Depends:` to `Imports:` to avoid unnecessary
  namespace clashes and encourage explicit fmesher declarations in user
  code (version `2.14.0.9003`). Temporarily keep re-exporting `fm_int`
  and `fm_pixels` to avoid breaking outdated packages, for CRAN bugfix
  release `2.14.1`.

## inlabru 2.14.0

CRAN release: 2026-03-08

### New features

- Add
  [`bru_names()`](https://inlabru-org.github.io/inlabru/reference/bru_names.md)
  method for extracting the inlabru standardised names of fixed effects,
  latent components, and hyperparameters from a fitted `bru` object.
  (version `2.13.0.9011`)
- Add
  [`bm_logitaverage()`](https://inlabru-org.github.io/inlabru/reference/bm_logitaverage.md)
  mapper, for weighted logit-averages (version `2.13.0.9010`)
- Add
  [`bm_reparam()`](https://inlabru-org.github.io/inlabru/reference/bm_reparam.md)
  mapper, for fixed-matrix reparameterisations of existing mappers
  (version `2.13.0.9019`)
- Add
  [`bru_comp_env_extra()`](https://inlabru-org.github.io/inlabru/reference/bru_comp_env_extra.md)
  getter and setter methods for the `env_extra` element of
  [`bru_comp()`](https://inlabru-org.github.io/inlabru/reference/bru_comp.md)
  objects, to allow storing extra information for use in component
  definitions. Mostly for internal use, but may be used by external
  packages needing to store and access special data. (version
  `2.13.0.9035`)
- Add
  [`bru_input_text()`](https://inlabru-org.github.io/inlabru/reference/bru_input_text.md)
  methods for extracting the component input information as text,
  e.g. for use in error messages (version `2.13.0.9036`)

### Updated features

- Regenerated `gorillas_sf` data set, with new mesh better suited for
  modelling, with a regular interior triangular mesh. (version
  `2.13.0.9034`)
- Regenerated `mexdolphin_sf` data set, with new mesh. (version
  `2.13.0.9035`)
- Add
  [`bru_get_mapper()`](https://inlabru-org.github.io/inlabru/reference/bru_get_mapper.md)
  support for `inla.cgeneric` objects, and update the support for
  `inla.rgeneric`, to standardise where to store a pre-constructed
  mapper, so that external packages will no longer need to add their own
  [`bru_get_mapper()`](https://inlabru-org.github.io/inlabru/reference/bru_get_mapper.md)
  methods for their r/cgeneric sub-classes, if the mapper is
  pre-computed. (version `2.13.0.9005`)
- Remove unused `group` and `replicate` parts of each latent component,
  speeding up component evaluation and linearisation. (version
  `2.13.0.9006`)
- Allow index extraction from duplicated tags in
  [`bru_index()`](https://inlabru-org.github.io/inlabru/reference/bru_index.md),
  returning all matching indices. (version `2.13.0.9014`)
- Code refactor and storage updates for
  [`bru_obs()`](https://inlabru-org.github.io/inlabru/reference/bru_obs.md)
  and
  [`bru_comp()`](https://inlabru-org.github.io/inlabru/reference/bru_comp.md)
  objects to prepare for future extensions. For 2.14, the compatibility
  option `bru_compat_pre_2_14_enable` is set to `TRUE` (version
  `2.13.0.9015`, `2.13.0.9017`, `2.13.0.9018`)
- Add `block_response` element to `bru_obs(aggregate_input)` list input,
  to allow specifying a block-wise response variable to match for
  aggregated predictors. (version `2.13.0.9016`)
- Deprecate the `bru_obs(allow_combine)` argument with its logical
  inverse, `bru_obs(is_rowwise)`, to make the meaning clearer. (version
  `2.13.0.9017`)
- Speed up `control.gcpo` friends structure generation for
  [`bru_obs()`](https://inlabru-org.github.io/inlabru/reference/bru_obs.md)
  with large numbers of samplers or observations. (version `2.13.0.9028`
  with further speedup in `2.13.0.9029`)
- Only enable `control.gcpo` computation in the final
  [`inla()`](https://rdrr.io/pkg/INLA/man/inla.html) iteration (version
  `2.13.0.9030`)
- Only enable `control.gcpo` for “cp” models if `control.gcpo` list
  input is supplied (version `2.13.0.9031`)
- Set option `bru_compress_cp = FALSE` by default, to avoid dense
  precision matrix blocks (version `2.13.0.9031`)
- Convert
  [`sample.lgcp()`](https://inlabru-org.github.io/inlabru/reference/sample.lgcp.md)
  to `sf` internals and output format (version `2.13.0.9033`).

### Bug fixes

- Fix logic for `bru_set_missing<inla.surv>()` to correctly handle the
  [`inla.surv()`](https://rdrr.io/pkg/INLA/man/surv.html) class
  variations (version `2.13.0.9002`)
- Handle changing number of hyperparameters in inla tracing (version
  `2.13.0.9004`)
- Robustify internal `extended_bind_rows()` method for unifying XY/XYZ
  `sf` coordinate columns (version `2.13.0.9009`)
- Robustify sf coordinate handling in
  [`generate()`](https://inlabru-org.github.io/inlabru/reference/generate.md)
  and
  [`predict()`](https://rspatial.github.io/terra/reference/predict.html)
  (version `2.13.0.9027`)
- Robustify sf coordinate handling in `bru_obs_family_cp()` for unifying
  XY/XYZ coordinate columns (version `2.13.0.9033`)
- Add logic for
  [`bru_fill_missing()`](https://inlabru-org.github.io/inlabru/reference/bru_fill_missing.md)
  to handle missing values in `sf` input data, so that
  e.g. `bru_fill_missing(input, input, input$values)` works as expected
  (version `2.13.0.9022`)
- Handle `bru` object upgrades from 2.12.0; the 2.12.0.9014 upgrade step
  needed to introduce the `"bru_obs"` class name earlier, and upgrades
  through 2.13.0.9017, that introduced a new predictor expression
  storage system (version `2.13.0.9024`)
- Fix bug in
  [`bru_obs_control_gcpo()`](https://inlabru-org.github.io/inlabru/reference/bru_obs_methods.md)
  for multi-observation models. (version `2.13.0.9025`)

## inlabru 2.13.0

CRAN release: 2025-07-09

### New features

- Allow the optional `weights` argument to
  `bru_obs(family = "cp")`/[`lgcp()`](https://inlabru-org.github.io/inlabru/reference/lgcp.md)
  to contain individual point observation weights, so that the `eta`
  contribution to the log-likelihood is `sum(weights * eta)` instead of
  `sum(eta)` (version `2.12.0.9006`)
- Automatically detect purely additive linear models (version
  `2.12.0.9014`)
- Add experimental predictor aggregation helper feature
  `bru_obs(..., aggregate = ..., aggregate_input = ...)` to simplify
  specification of models with aggregation as the final step of the
  predictor evaluation (version `2.12.0.9013`, bugfix in `2.12.0.9016`).
  Includes support for constructing the aggregation information, via
  `domain`,`samplers`, or precomputed `ips` (version `2.12.0.9022`)
- Add
  [`bru_set_missing()`](https://inlabru-org.github.io/inlabru/reference/bru_set_missing.md)
  method for setting missing values in the `bru_obs` data, e.g. for use
  in cross-validation or prior sampling (version `2.12.0.9024`)

### Updates features

- Add automated support for INLA models with hidden states, beyond the
  “bym” and “bym2” models (version `2.12.0.9002`)
- Add automatic mapper support for INLA lattice models “rw2d”,
  “rw2diid”, and “matern2d” (version `2.12.0.9004`)
- Add `inputs` data to the `bru_info` object, so that the component
  inputs can be pre-evaluated before automated mapper construction, and
  avoiding duplicate evaluation in
  [`iinla()`](https://inlabru-org.github.io/inlabru/reference/iinla.md).
  This can halve the pre-processing time for large spatial and
  spatio-temporal models (version `2.12.0.9010`)
- Add
  [`bru_log()`](https://inlabru-org.github.io/inlabru/reference/bru_log.md)
  data for warnings and errors reported by inlabru (version
  `2.12.0.9011`)
- Rename `bru_like` object class to `bru_obs` (version `2.12.0.9017`)
  and temporarily re-reintroduce
  [`like_list()`](https://inlabru-org.github.io/inlabru/reference/bru_obs.md)
  and
  [`bru_like_list()`](https://inlabru-org.github.io/inlabru/reference/bru_obs.md)
  as aliases for
  [`bru_obs_list()`](https://inlabru-org.github.io/inlabru/reference/bru_obs.md)
  (version `2.12.0.9020`)
- Add `quantile` argument to
  [`spde.posterior()`](https://inlabru-org.github.io/inlabru/reference/spde.posterior.md)
  for controlling the credible interval calculations, like
  [`materncov.bands()`](https://inlabru-org.github.io/inlabru/reference/materncov.bands.md),
  which is now also exported (version `2.12.0.9019`)
- Remove `dic` and `waic` from default `control.compute` `bru_options`,
  to match INLA defaults (version `2.12.0.9024`)

### New and updated mapper features

- The standard mapper class names have been shortened from
  `bru_mapper_<type>` to `bm_<type>`, to make code using them more
  readable. Objects of the old class names will be converted to the new
  classes internally, so that old stored objects will still work.
  Constructors of the form `bru_mapper_<type>()` call the corresponding
  new constructor `bm_<type>()` (version `2.12.0.9021`)
- Add
  [`bm_sum()`](https://inlabru-org.github.io/inlabru/reference/bm_sum.md)
  mapper, for automated adding the output of multiple mappers,
  optionally with a single common input (version `2.12.0.9001`)
- Add `interleaved` option to
  [`bm_repeat()`](https://inlabru-org.github.io/inlabru/reference/bm_repeat.md)
  to allow interleaved states for summation of a repeated mapper
  (version `2.12.0.9001`)
- Allow `n_block` in `input` argument to `bm_aggregate` and
  `bm_logsumexp` evaluation methods, overriding the optional mapper
  object setting (version `2.12.0.9005`)
- (Note: due to the difficulty of ensuring correct output ordering, and
  `fmesher` will refuse `character` block input from version `0.5.0`.
  Was: Allow `character` block information in `bm_aggregate` and
  `bm_logsumexp` mappers, from `fmesher` version `0.2.0.9017` (version
  `2.12.0.9013`))
- Expanded auto-detection of component sizes by checking `Cmatrix` and
  `graph` arguments, if present and `n` is `NULL` (version
  `2.12.0.9012`)
- Code refactor to expand `bm_list` mapper list handling, removing
  unnecessary method layers for component linearisation and
  simplification (version `2.12.0.9018`)

### Bugfixes and deprecations

- Give deprecation warnings for `Spatial` object inputs to
  [`bru_obs()`](https://inlabru-org.github.io/inlabru/reference/bru_obs.md),
  as maintaining the fallback support code is becoming increasingly
  time-consuming, and the `sf` package is now the recommended spatial
  data handling package. (version `2.12.0.9023`)
- Bugfix for
  [`generate.bru()`](https://inlabru-org.github.io/inlabru/reference/generate.md)
  for evaluation of expressions in the absence of `newdata` (version
  `2.12.0.9007`)
- Use `mid` locations for
  [`ibm_values()`](https://inlabru-org.github.io/inlabru/reference/ibm_values.md)
  for non-indexed `fm_mesh_1d` mapper (version `2.12.0.9009`)
- Deprecate the `include` and `exclude` arguments to
  [`predict()`](https://rspatial.github.io/terra/reference/predict.html)
  and
  [`generate()`](https://inlabru-org.github.io/inlabru/reference/generate.md)
  (version `2.12.0.9003`)
- Assumed that [`NROW()`](https://rdrr.io/r/base/nrow.html) on the main
  component gave the correct size for `group` and `replicate` in the
  component `_eval()` feature. Now uses
  [`ibm_n_output()`](https://inlabru-org.github.io/inlabru/reference/ibm_n_output.md)
  instead. Fixes
  [\#271](https://github.com/inlabru-org/inlabru/issues/271) (version
  `2.12.0.9015`)

## inlabru 2.12.0

CRAN release: 2024-11-21

### General changes

- Introduce
  [`bru_obs()`](https://inlabru-org.github.io/inlabru/reference/bru_obs.md)
  as a replacement to
  [`like()`](https://inlabru-org.github.io/inlabru/reference/bru_obs.md),
  to help avoid namespace clashes with
  e.g. [`data.table::like()`](https://rdrr.io/pkg/data.table/man/like.html)
  (version `2.11.1.9026`)
- Change logic for
  `bru_obs(allow_combine)`/[`like()`](https://inlabru-org.github.io/inlabru/reference/bru_obs.md)
  to allow user override, and add warnings for ambiguous cases (version
  `2.11.1.9011`)
- Add
  [`bru_response_size()`](https://inlabru-org.github.io/inlabru/reference/bru_response_size.md)
  method for extracting the response size for each observation `bru_obs`
  object (version `2.11.1.9013`)
- Add `sf` output format support for `sline` and `spoly` (version
  `2.11.1.9006`)
- Add `[` and `]` to special character set in
  [`bru_standardise_names()`](https://inlabru-org.github.io/inlabru/reference/bru_standardise_names.md)
  (version `2.11.1.9012`)
- Add
  [`bru_index()`](https://inlabru-org.github.io/inlabru/reference/bru_index.md)
  method for accessing predictor index information for sub-models, and a
  `tag` argument for
  [`bru_obs()`](https://inlabru-org.github.io/inlabru/reference/bru_obs.md)
  to identify individual sub-models by name, which is also propagated to
  lists of `bru_like` objects (version `2.11.1.9017` and `2.12.0`)
- Allow
  [`bru_mapper_multi()`](https://inlabru-org.github.io/inlabru/reference/bm_multi.md)
  sub-mappers to have non-zero offsets, that are added to generate the
  combined offset (version `2.11.1.9019`)
- Add general
  [`bru_mapper_fmesher()`](https://inlabru-org.github.io/inlabru/reference/bm_fmesher.md)
  mapper, for indexed mapping of all objects supporting
  [`fm_dof()`](https://inlabru-org.github.io/fmesher/reference/fm_dof.html)
  and
  [`fm_basis()`](https://inlabru-org.github.io/fmesher/reference/fm_basis.html)
  (version `2.11.1.9021`)
- Add
  [`bru_mapper_repeat()`](https://inlabru-org.github.io/inlabru/reference/bm_repeat.md)
  mapper, for automated single-mapper sums (version `2.11.1.9022`)

### Data set updates

- Convert `mrsea` to `sf` format (version `2.11.1.9009`)
- Convert the `shrimp` data set to `sf` format (version `2.11.1.9007`)
- Remove `seals_sp` data set due to excessive size (version
  `2.11.1.9010`)
- Replace the `mexdolphin` dataset with a function
  [`mexdolphin_sp()`](https://inlabru-org.github.io/inlabru/reference/mexdolphin_sf.md)
  to avoid `sp` data objects in the package (version `2.11.1.9015`)
- Replace the `gorillas` dataset with a function
  [`gorillas_sp()`](https://inlabru-org.github.io/inlabru/reference/gorillas_sf.md)
  to avoid `sp` data objects in the package (version `2.11.1.9016`)

### Namespace changes

- Move `sp` from `Imports` to `Suggests`. Component definitions using
  `coordinates` as input require either
  [`sp::coordinates`](https://edzer.github.io/sp/reference/coordinates.html)
  or `sp` having been already loaded with
  e.g. [`library(sp)`](https://github.com/edzer/sp/) (version
  `2.11.1.9003`)
- Remove `ggmap` support (version `2.11.1.9002`)
- Remove unnecessary `INLA` namespace loading in `ggplot` methods
  (version `2.11.1.9008`)
- Move `terra` from `Imports` to `Suggests` (version `2.11.1.9014`)
- Stop re-exporting `fmesher` methods (version `2.11.1.9020`)

### Deprecated methods

- Deprecated (since 2.8.0) method `is.inside()` have been removed. Use
  [`fmesher::fm_is_within()`](https://inlabru-org.github.io/fmesher/reference/fm_is_within.html)
  instead.
- Deprecated (since 2.7.0) `bru_mapper.default()` to define new mapper
  classes has been removed. Use
  [`bru_mapper_define()`](https://inlabru-org.github.io/inlabru/reference/bru_mapper.md)
  instead.
- Deprecated (since 2.6.0) `bru_mapper_offset()` method has been
  removed. Use
  [`bru_mapper_const()`](https://inlabru-org.github.io/inlabru/reference/bm_const.md)
  instead.

### Internal changes

- Remove unneeded `"list"` class inheritance from solitary classes
  (version `2.11.1.9001`)
- Expand the `summary` and `print` method class coverage (version
  `2.11.1.9002`)
- Reduced the amount of diagnostic messages in
  [`bru_safe_inla()`](https://inlabru-org.github.io/inlabru/reference/bru_safe_inla.md)
  (version `2.11.1.9005`)
- Change name of `component` and `component_list` methods to
  `bru_component` and `bru_component_list` (version `2.11.1.9026`)

## inlabru 2.11.1

CRAN release: 2024-07-01

### Bug fixes

- Fix documentation link issues not spotted by local and github package
  checks
- Work around invalid geometry in map data used by the Spatially Varying
  Coefficients vignette

## inlabru 2.11.0

### New features

- Add support for the `scale` parameter to
  [`like()`](https://inlabru-org.github.io/inlabru/reference/bru_obs.md).
- Add support for special response objects like
  [`inla.mdata()`](https://rdrr.io/pkg/INLA/man/mdata.html) and
  [`inla.surv()`](https://rdrr.io/pkg/INLA/man/surv.html), when INLA
  version `> 24.06.02` (for `mdata`) or `> 24.06.26` (for `surv`) are
  available (version `2.10.1.9011`)
- New `toypoints` example data set, for basic modelling examples
  (version `2.10.1.9003`)

### Feature updates

- Updated convergence plots, reducing random effect aspects to summary
  statistics, improving speed and visual coherence (version
  `2.10.1.9004`)
- Add options to
  [`bru_convergence_plot()`](https://inlabru-org.github.io/inlabru/reference/bru_convergence_plot.md)
  to control the number of iterations shown, and optionally show the
  initial values that are stored from this version (version
  `2.10.1.9005`)
- Switch timing mechanism from
  [`Sys.time()`](https://rdrr.io/r/base/Sys.time.html) to
  [`proc.time()`](https://rdrr.io/r/base/proc.time.html) to capture CPU
  time instead of elapsed clock time. Added
  [`bru_timings()`](https://inlabru-org.github.io/inlabru/reference/bru_timings.md)
  method to extract the timings safely from a fitted `bru` object
  (version `2.10.1.9007` and `2.10.1.9010`)
- Add verbosity level information to the bru log data structure,
  allowing filtered log extraction and more flexible log display
  (version `2.10.1.9012`)

### Bug fixes

- Fix regression bug in `"bym"` model support, where the latent state
  size wasn’t correctly handled by the mapper system (version
  `2.10.1.9002`)
- Add filter to limit mapper construction to only the components used in
  the predictor expression, to avoid unused components breaking the
  initialisation. This allows easier testing of multi-likelihood models
  (version `2.10.1.9006`)
- Improved backwards compatibility support for `sp` data input for
  `family = "cp"` (version `2.10.1.9008`)

### Deprecated methods

- Deprecated (since 2.9.0) method `ipoints(samplers, domain)` is no
  longer available. Use `fmesher::fm_int(domain, samplers)` instead.
- The `allow_latent`, `include_latent` arguments to
  [`like()`](https://inlabru-org.github.io/inlabru/reference/bru_obs.md)
  have been deprecated in favour of the general
  [`bru_used()`](https://inlabru-org.github.io/inlabru/reference/bru_used.md)
  framework, that auto-detects what component effects and latent effects
  are used by a predictor expression.
- The deprecated (since 2.8.0) `cprod()` method now gives a warning and
  will be removed in a future version. Use
  [`fmesher::fm_cprod()`](https://inlabru-org.github.io/fmesher/reference/fm_cprod.html)
  instead.
- The `integration_weight_aggregation` method has been removed
  (deprecated since 2.8.0). Use
  [`fmesher::fm_vertex_projection()`](https://inlabru-org.github.io/fmesher/reference/fm_vertex_projection.html)
  instead.
- The `mesh_triangle_integration` method has been removed (deprecated
  since 2.8.0). Use
  [`fmesher::fm_int()`](https://inlabru-org.github.io/fmesher/reference/fm_int.html)
  instead.
- Old use of `bru_mapper.default()` to define new mapper classes has
  been disabled (deprecated since 2.7.0). Use
  [`bru_mapper_define()`](https://inlabru-org.github.io/inlabru/reference/bru_mapper.md)
  instead.
- Deprecated (since 2.8.0) methods `is.inside()`,
  `vertices.inla.mesh()`, and `pixels()` have been disabled. Use
  [`fmesher::fm_is_within()`](https://inlabru-org.github.io/fmesher/reference/fm_is_within.html),
  [`fmesher::fm_vertices()`](https://inlabru-org.github.io/fmesher/reference/fm_vertices.html),
  and
  [`fmesher::fm_pixels()`](https://inlabru-org.github.io/fmesher/reference/fm_pixels.html)
  instead.

## inlabru 2.10.1

CRAN release: 2023-12-21

### Feature updates

- Add new web article on ZIP/ZAP models (zero inflated Poisson models)
  where the non-zero probability is modelled with a separate predictor
  from the predictor for the Poisson model parameter:
  <https://inlabru-org.github.io/inlabru/articles/zip_zap_models.html>.
  (Thanks to Dmytro Perepolkin)

### Bug fixes and dependency simplification

- Remove dependence on the `ggpolypath` package, and the
  `ggplot2::fortify.SpatialPolygons/DataFrame()` methods that were
  deprecated in `ggplot2` version `3.4.4`. Code using
  [`gg.SpatialPolygons()`](https://inlabru-org.github.io/inlabru/reference/gg.Spatial.md)
  together with
  [`coord_fixed()`](https://ggplot2.tidyverse.org/reference/coord_fixed.html)/[`coord_equal()`](https://ggplot2.tidyverse.org/reference/coord_fixed.html)
  for coordinate axis control needs to use
  [`coord_sf()`](https://ggplot2.tidyverse.org/reference/ggsf.html)
  instead.
- Detect the need for vectorised parameters in
  `bru_forward_transformation` to allow `bru_mapper_marginal` to be
  applied with e.g. spatially varying parameters. (version
  `2.10.0.9001`)
- Detect `terra` version `>= 1.7-66` that removes the need for detecting
  special cases (`nrow(where) == 1` and `terra::nlyr(data) == 1`).
  Workaround code used for versions `< 1.7-66`. (version `2.10.0.9002`)
  (Thanks to Robert J. Hijmans)

## inlabru 2.10.0

CRAN release: 2023-10-29

### Feature updates

- Add new
  [`ibm_simplify()`](https://inlabru-org.github.io/inlabru/reference/ibm_simplify.md)
  generic to handle mapper simplification more generally; needed to
  properly support non-linear component mappers. (version `2.9.0.9004`)
- Add new
  [`bru_mapper_marginal()`](https://inlabru-org.github.io/inlabru/reference/bm_marginal.md)
  mapper class that can be used as part of component mapper pipelines.
  (version `2.9.0.9004`)
- Add new
  [`ibm_eval2()`](https://inlabru-org.github.io/inlabru/reference/ibm_eval2.md)
  generic that computes both evaluation and Jacobian, avoiding
  double-computing of the Jacobian, when practical. (version
  `2.9.0.9005`)
- Add new
  [`bru_timings_plot()`](https://inlabru-org.github.io/inlabru/reference/bru_timings_plot.md)
  function that plots the time used for each nonlinear iteration
  (version `2.9.0.9007`)
- Speed up
  [`bru_fill_missing()`](https://inlabru-org.github.io/inlabru/reference/bru_fill_missing.md)
  (by orders of magnitude) by changing method for finding the nearest
  available data point. (version `2.9.0.9011`)
- Add new
  [`bru_mapper_shift()`](https://inlabru-org.github.io/inlabru/reference/bm_shift.md)
  mapper class that works like
  [`bru_mapper_scale()`](https://inlabru-org.github.io/inlabru/reference/bm_scale.md)
  but for additive shifts instead of multiplicative scaling. (version
  `2.9.0.9012`)
- Added more checks for invalid component or predictor evaluations, to
  help catch user errors sooner, and with more informative messages.
  (version `2.9.0.9013`)
- Expand `bru_mapper_matrix`, previously used only for component
  `model = "fixed"`, to allow integer indexing in addition to the
  previous factor/character-only indexing. (version `2.9.0.9014`)

### Bug fixes

- The `is_linear` flag wasn’t correctly set for `bru_mapper_logsumexp`
  mappers. Since previous versions did not accept non-linear component
  mappers, this is unlikely to have affected any user code. (Fixed in
  version `2.9.0.9001`)
- Improved error messages for missing or incomplete LGCP domain
  specification. (version `2.9.0.9002` and `2.9.0.9006`)
- Allow `NULL` in automatic component usage detection. (version
  `2.9.0.9003`)
- Corrected the crs information for `gorillas$plotsample$counts` and
  `gorillas_sf$plotsample$counts` from `+units=m` to `+units=km`.
  (version `2.9.0.9010`) The geometry information in `counts` is
  unlikely to have been used in examples or analysis code, as the
  problem would have been immediately obvious; plotting or other
  geometric operations that use the crs information would heve been
  completely wrong, and is only detected now that more code uses the crs
  information at all. Thanks to Dmytro Perepolkin for reporting in issue
  [\#205](https://github.com/inlabru-org/inlabru/issues/205)
- Fix problem in
  [`bru_fill_missing()`](https://inlabru-org.github.io/inlabru/reference/bru_fill_missing.md)
  for cases where the input data object also has missing values.
  (version `2.9.0.9011`)
- Make
  [`eval_spatial()`](https://inlabru-org.github.io/inlabru/reference/eval_spatial.md)
  transform the `where` coordinates to the same crs as the input data,
  for `SpatRaster` and `sf` inputs, to allow different crs
  specifications. (version `2.9.0.9012`)

## inlabru 2.9.0

CRAN release: 2023-08-28

### Feature updates

- Conversion of code to use `fmesher` for mesh and geometry handling;
  the interface supports existing objects and methods. See
  <https://inlabru-org.github.io/fmesher/articles/inla_conversion.html>
  for more information.

- General speed improvements, see below for details.

- Added
  [`gg.sf()`](https://inlabru-org.github.io/inlabru/reference/gg.sf.md)
  method.

- Add experimental support for `stars` via
  [`eval_spatial()`](https://inlabru-org.github.io/inlabru/reference/eval_spatial.md).
  (version `2.8.0.9007`)

- Move the `sp` package from ‘Depends’ to ‘Imports’. This means that
  user code should either use `sp::` or
  [`library("sp")`](https://github.com/edzer/sp/) to access `sp`
  methods. The
  [`bru_safe_sp()`](https://inlabru-org.github.io/inlabru/reference/bru_safe_sp.md)
  helper function can be used to check for a safe `sp` package
  configuration during the transition from `rgdal` to `sf`, and is only
  needed if you may run on systems with `sp` installations older than
  “2.0-0” or with `sp::get_evolution_status() < 2`. (version `2.8.2011`)

- Now preserves the previous log output when using
  [`bru_rerun()`](https://inlabru-org.github.io/inlabru/reference/bru.md),
  and
  [`bru_log()`](https://inlabru-org.github.io/inlabru/reference/bru_log.md)
  is now a set of S3 methods, supporting extracting the full inlabru log
  as well `bru`-object specific logs (version `2.8.0.9008`).

  Note: From version `2.9.0`, use
  [`bru_log()`](https://inlabru-org.github.io/inlabru/reference/bru_log.md)
  to access the global log, and `bru_log(fit)` to access a stored
  estimation log.

  Up to version `2.8.0`,
  [`bru_log()`](https://inlabru-org.github.io/inlabru/reference/bru_log.md)
  was a deprecated alias for
  [`bru_log_message()`](https://inlabru-org.github.io/inlabru/reference/bru_log_message.md).
  When running on `2.8.0` or earlier, use `bru_log_get()` to access the
  global log, and `cat(fit$bru_iinla$log, sep = "\n")` to print a stored
  estimation object log.

### Bug fixes and speed improvements

- Covariate object component inputs of type `SpatialPolygonsDataFrame`
  were not automatically passed on to
  [`eval_spatial()`](https://inlabru-org.github.io/inlabru/reference/eval_spatial.md).
  The logic has now changed so that any object with a
  [`eval_spatial()`](https://inlabru-org.github.io/inlabru/reference/eval_spatial.md)
  method will trigger a call to
  [`eval_spatial()`](https://inlabru-org.github.io/inlabru/reference/eval_spatial.md).
  See
  [`?input_eval`](https://inlabru-org.github.io/inlabru/reference/inlabru-deprecated.md)
  for further information. (version `2.8.0.9001`)
- [`fm_crs_is_null()`](https://inlabru-org.github.io/fmesher/reference/fm_crs_is_null.html),
  [`fm_transform()`](https://inlabru-org.github.io/fmesher/reference/fm_transform.html)
  now supports oblique `fm_crs` CRS objects, and
  [`is.na()`](https://rdrr.io/r/base/NA.html) methods for the `fm_crs`
  and `inla.CRS` classes have been added. (version `2.8.0.9003`)
- Significant speed up
  [`predict()`](https://rspatial.github.io/terra/reference/predict.html)
  by using `quantile(..., names = FALSE)`. (version `2.8.0.9004`)
- Improved
  [`row_kron()`](https://inlabru-org.github.io/inlabru/reference/inlabru-deprecated.md)
  code, causing speedups of a factor 2-30 in randomised test cases.
  (version `2.8.0.9005`)
- Removed incorrect code for `sf` method for
  [`eval_spatial()`](https://inlabru-org.github.io/inlabru/reference/eval_spatial.md),
  causing failure when extracting from multiple layers in a single call.
  (version `2.8.0.9007`)
- Improved handling of posterior sample variable extraction in
  [`generate()`](https://inlabru-org.github.io/inlabru/reference/generate.md)
  and
  [`predict()`](https://rspatial.github.io/terra/reference/predict.html).
  Now much faster for large models. (version `2.8.0.9009`)
- Fixed linearisation issue when using only the `*_latent` form of a
  component. (version `2.8.0.9015`)
- Workaround for equivalent but textually different CRS/WKT information
  in
  [`bru_fill_missing()`](https://inlabru-org.github.io/inlabru/reference/bru_fill_missing.md).
  (version `2.8.0.9016`, fixes
  [\#200](https://github.com/inlabru-org/inlabru/issues/200))

### Deprecation of old functions

- `eval_SpatialDF` removed, deprecated since `2.8.0`. See `eval_spatial`
  instead.
- `stransform`, `ibm_amatrix`, `ibm_valid_input` removed, deprecated
  since `2.7.0`. See `fm_transform` and `ibm_jacobian` instead.
- `bru_mapper_offset`, deprecated since `2.6.0` now returns a pure
  `bru_mapper_const` object, and all `bru_mapper_offset` `ibm_*` methods
  have been removed.
- `init.tutorial` removed, deprecated since `2.5.0`
- `generate.inla` and `predict.inla` removed, deprecated since `2.1.0`

## inlabru 2.8.0

CRAN release: 2023-06-20

### Feature updates

- The iterative inla method has been given both sharper internal
  [`inla()`](https://rdrr.io/pkg/INLA/man/inla.html) optimisation
  criteria for the iterations (thanks to Haavard Rue), *and* a more
  relaxed nonlinear iteration stopping criterion; the default
  `bru_method$rel_tol` values has been changed from 1 to 10 percent
  change. The iterations are terminated when all latent and
  hyper-parameter mode changes fullfil `|change|/SD < rel_tol`, and the
  non-linear line search is inactive. This seems to strike a useful
  balance between the different optimisation criteria, allowing the
  iterations to converge faster and also detect that convergence sooner.

- The logic for which components are needed for a predictor expression
  (in
  [`like()`](https://inlabru-org.github.io/inlabru/reference/bru_obs.md)
  or
  [`generate()`](https://inlabru-org.github.io/inlabru/reference/generate.md)/[`predict()`](https://rspatial.github.io/terra/reference/predict.html))
  has been updated to when possible extract the list of components from
  the expression itself. The user can override this default if
  necessary, using the `include`/`exclude` arguments.

  The
  [`bru_used()`](https://inlabru-org.github.io/inlabru/reference/bru_used.md)
  methods are used to guess the needed component names, applied to the
  right-hand side of the `formula` arguments. The `allow_latent`
  argument to
  [`like()`](https://inlabru-org.github.io/inlabru/reference/bru_obs.md)
  has been deprecated in favour of `include_latent` (by default
  auto-detected for use of `_latent` and `_eval`).

  The internal information storage is handled by the new
  [`bru_used()`](https://inlabru-org.github.io/inlabru/reference/bru_used.md)
  methods, that can also be used directly by the user and supplied via
  the `used` argument to
  [`like()`](https://inlabru-org.github.io/inlabru/reference/bru_obs.md)/[`generate()`](https://inlabru-org.github.io/inlabru/reference/generate.md)/[`predict()`](https://rspatial.github.io/terra/reference/predict.html).

- Add
  [`fm_int()`](https://inlabru-org.github.io/fmesher/reference/fm_int.html)
  integration methods, replacing the old `ipmaker()` and `ipoints()`
  methods. Supports both `sf` and `sp` sampler objects.

- Add
  [`fm_pixels()`](https://inlabru-org.github.io/fmesher/reference/fm_pixels.html)
  methods for gridded points. The old `pixels()` method now calls
  `fm_pixels(..., format = "sp")`

- `eval_spatial` support for sf objects (for point-in-polygon data
  lookups)

- Allow precomputed spatial covariates in the data for point process
  observations

- Add `edge|int|ext.linewidth` arguments to `gg.inla.mesh`
  [\#188](https://github.com/inlabru-org/inlabru/issues/188)

- Rename the
  [`predict()`](https://rspatial.github.io/terra/reference/predict.html)
  and
  [`generate()`](https://inlabru-org.github.io/inlabru/reference/generate.md)
  `data` arguments to `newdata`, for better compatibility with other
  [`predict()`](https://rspatial.github.io/terra/reference/predict.html)
  methods. The old argument name will still be accepted, but give a
  warning. Code that does not name the `data` argument is not affected.

- Note: Coordinate names for `Spatial*` objects have been inconsistently
  available in the predictor expression evaluation. However, due to how
  internal conversions might inadvertently change these names, they can
  not be relied on, and they are no longer being made available to the
  predictor expression. As a side effect, this change also speeds up
  some [`bru()`](https://inlabru-org.github.io/inlabru/reference/bru.md)
  runs by around a factor 2, since it avoids converting the `Spatial*`
  to a regular `data.frame` in time-sensitive core evaluation code.

  If you need access to the raw coordinate values, use explicit calls to
  `sp::coordinates(.data.)` (e.g. for custom spatial covariate
  evaluation.). When possible, use the built-in covariate evaluation
  method,
  [`eval_spatial()`](https://inlabru-org.github.io/inlabru/reference/eval_spatial.md),
  either implicitly with `comp(covariate, ...)` or explicitly,
  `comp(eval_spatial(covariate, where = .data.), ...)`, that handles
  `crs` information correctly. Also consider transitioning from `sp` to
  `sf` data storage, using `geometry` instead of raw coordinates.

### Bug and dependency updates

- Remove `rgdal` and `maptools` dependencies
  [\#178](https://github.com/inlabru-org/inlabru/issues/178)
- Add
  [`bru_safe_sp()`](https://inlabru-org.github.io/inlabru/reference/bru_safe_sp.md)
  to check if `sp` can be used safely (checks `rgdal` availability and
  `sp` evolution status, optionally forcing use of `sf`)
  [\#178](https://github.com/inlabru-org/inlabru/issues/178)
- Remove PROJ4 support
  [\#178](https://github.com/inlabru-org/inlabru/issues/178)
- Change `rgl.*` functions to `*3d`. Thanks to Duncan Murdoch
  [\#181](https://github.com/inlabru-org/inlabru/issues/181)
- Speed up `ibm_jacobian.bru_mapper_harmonics` for large models
- Add workarounds for inconsistent polygon orientation resulting from
  `sf::st_*` calls that don’t account for the `geos` canonical
  representation being CW, whereas the canonical Simple Features
  representation being CCW. See
  <https://github.com/r-spatial/sf/issues/2096>

## inlabru 2.7.0

CRAN release: 2022-12-02

### Feature overview

- Added support for `sf` and `terra` inputs to most methods
- Expanded geometry and mesh handling methods
- Expanded
  [`bru_mapper()`](https://inlabru-org.github.io/inlabru/reference/bru_mapper.md)
  system
- Added convergence diagnostics plot with
  [`bru_convergence_plot()`](https://inlabru-org.github.io/inlabru/reference/bru_convergence_plot.md)

### Feature details

- Allow `NA` input for default 1D mappers to generate effect zero, like
  in [`inla()`](https://rdrr.io/pkg/INLA/man/inla.html).

- New and expanded methods
  [`fm_crs()`](https://inlabru-org.github.io/fmesher/reference/fm_crs.html),
  [`fm_CRS()`](https://inlabru-org.github.io/fmesher/reference/fm_CRS_sp.html),
  [`fm_transform()`](https://inlabru-org.github.io/fmesher/reference/fm_transform.html),
  [`fm_ellipsoid_radius()`](https://inlabru-org.github.io/fmesher/reference/fm_crs_wkt.html),
  and
  [`fm_length_unit()`](https://inlabru-org.github.io/fmesher/reference/fm_crs_wkt.html)
  to further support `sf` objects. The
  [`fm_crs()`](https://inlabru-org.github.io/fmesher/reference/fm_crs.html)
  extraction method also supports `terra` objects.

- [`bru_fill_missing()`](https://inlabru-org.github.io/inlabru/reference/bru_fill_missing.md)
  now supports `terra` `SpatRaster` data and and `sf` locations.

- New experimental methods
  [`fm_evaluator()`](https://inlabru-org.github.io/fmesher/reference/fm_evaluate.html)
  and
  [`fm_evaluate()`](https://inlabru-org.github.io/fmesher/reference/fm_evaluate.html),
  replacing the `INLA` `inla.mesh.projector` and `inla.mesh.project`
  methods.

- Experimental integration support for sphere and globe meshes.

- Allow `sf` input to `family="cp"` models.

- Further
  [`bru_mapper()`](https://inlabru-org.github.io/inlabru/reference/bru_mapper.md)
  method updates;

  - Deprecated `ibm_amatrix()` and
    [`names()`](https://rspatial.github.io/terra/reference/names.html)
    methods, replaced by
    [`ibm_jacobian()`](https://inlabru-org.github.io/inlabru/reference/ibm_jacobian.md)
    and
    [`ibm_names()`](https://inlabru-org.github.io/inlabru/reference/ibm_names.md).
  - Introduced
    [`bru_mapper_pipe()`](https://inlabru-org.github.io/inlabru/reference/bm_pipe.md),
    used to link mappers in sequence.
  - Introduced
    [`bru_mapper_aggregate()`](https://inlabru-org.github.io/inlabru/reference/bm_aggregate.md)
    and
    [`bru_mapper_logsumexp()`](https://inlabru-org.github.io/inlabru/reference/bm_logsumexp.md),
    used for blockwise weighted sums and log-sum-exp mappings,
    `output[k] = sum(weights[block==k]*state[block==k])))` and
    `output[k] = log(sum(weights[block==k]*exp(state[block==k])))`, with
    optional weight normalisation within each block. Allows providing
    the weights as log-weights, and uses block-wise shifts to avoid
    potential overflow.
  - `summary` methods for `bru_mapper` objects
    ([`summary.bru_mapper()`](https://inlabru-org.github.io/inlabru/reference/bm_summary.md))
  - Removed `methods` argument from
    [`bru_mapper_define()`](https://inlabru-org.github.io/inlabru/reference/bru_mapper.md).
    Implementations should register S3 methods instead.

### Bug fixes

- Remove unused `spatstat.core` dependency. Fixes
  [\#165](https://github.com/inlabru-org/inlabru/issues/165)
- Fixed issue with plain mapper evaluation in the
  [`ibm_eval.default()`](https://inlabru-org.github.io/inlabru/reference/ibm_eval.md)
  and `ibm_eval.bru_mapper_collect()` methods, where they would return
  zeros instead of the intended values. The main component evaluation
  and estimation code was not directly affected as that is based on the
  [`bru_mapper_multi()`](https://inlabru-org.github.io/inlabru/reference/bm_multi.md)
  class methods that rely on the Jacobians instead. The bug would
  therefore mainly have impacted the future, not yet supported nonlinear
  mapper extensions.
- Fix for `eval_spatial.SpatRaster`; Work around inconsistent logic in
  `terra::extract(..., layer)` when `length(layer)==1` or
  `nrow(where)==1`. Fixes
  [\#169](https://github.com/inlabru-org/inlabru/issues/169)
- Add `indexed` logical option to
  [`bru_mapper_factor()`](https://inlabru-org.github.io/inlabru/reference/bm_factor.md),
  to allow factor inputs to be mapped to index values, as needed for
  `group` and `replicate`. Fixes
  [\#174](https://github.com/inlabru-org/inlabru/issues/174)

## inlabru 2.6.0

CRAN release: 2022-10-24

### Features

- Add `bru_get_mapper` generic, and associated methods for `inla.spde`
  and `inla.rgeneric` objects. This allows `inlabru` to automatically
  extract the appropriate `bru_mapper` object for each model component,
  and can be used as a hook by external packages implementing new INLA
  object classes.

- Add a `weights` argument for
  [`like()`](https://inlabru-org.github.io/inlabru/reference/bru_obs.md),
  for likelihood-specific log-likelihood weights, passed on to the
  [`INLA::inla()`](https://rdrr.io/pkg/INLA/man/inla.html) weights
  argument. Evaluated in the data context.

- The `<component>_eval()` methods available in predictor expressions
  now handle optional scaling weights, like in ordinary component effect
  evaluation.

- Add `terra` support for covariate inputs

- The component `*_layer` arguments are now evaluated in the data
  context, to allow dynamic layer selection for spatial raster
  covariates. A new generic
  [`eval_spatial()`](https://inlabru-org.github.io/inlabru/reference/eval_spatial.md)
  provides support for grid/pixel based `Spatial*DataFrame` evaluation,
  and `SpatRaster`. Expanded support is in progress.

- New vignettes on the `bru_mapper` system, `component` definitions, and
  `prediction_scores`

- General overhaul of the `bru_mapper` and linearised predictor system,
  to prepare for new features.

  - Add `ibm_eval` generic for evaluating mappers for given states.
  - Add `bru_mapper_taylor`, used as an internal mapper for linearised
    mappers. This and `ibm_eval` is aimed at future support for
    nonlinear mappers. Associated new generic methods:
    `ibm_{is_linear,jacobian,linear}`.
  - New mapper implementations should use `ibm_jacobian` instead of
    `ibm_amatrix`. This allows defining a linearised mapper via
    `ibm_eval(input, state0) + ibm_jacobian(input, state0) %*% (state - state0)`.
  - New mapper class `bru_mapper_const`, which replaces
    `bru_mapper_offset`. `bru_mapper_offset` is now deprecated and will
    produce warnings.

### Bug fixes

- Enable both datum/ensemble container for ellipsoid information, to
  support `epsg:4326`. Fixes
  [\#154](https://github.com/inlabru-org/inlabru/issues/154)
- Make duplicated component names an error instead of a warning. Relates
  to [\#155](https://github.com/inlabru-org/inlabru/issues/155)
- Fix `Tsparse` assumptions in `row_kron` to prepare for Matrix `1.5-2`.
  Fixes [\#162](https://github.com/inlabru-org/inlabru/issues/162)

## inlabru 2.5.3

CRAN release: 2022-09-05

### Features

- Add `bru_mapper_harmonics` mapper for `cos` and `sin` basis sets.
- Allow
  [`predict()`](https://rspatial.github.io/terra/reference/predict.html)
  input data to be be a list.
- Allow arbitrary quantile summaries in
  [`predict()`](https://rspatial.github.io/terra/reference/predict.html)
- Remove `cv`, `var`, `smin`, `smax` summaries from
  [`predict()`](https://rspatial.github.io/terra/reference/predict.html)
- Add `mean.mc_std_err` and `sd.mc_std_err` output to
  [`predict()`](https://rspatial.github.io/terra/reference/predict.html)
- Add `robins_subset` data set and associated variable coefficient web
  vignette

### Bug fixes

- Propagate multi-likelihood A-matrix information instead of
  recomputing. Fixes iteration issue for bym2 and other
  `bru_mapper_collect` models.
- Turn on predictor summaries during iterations to allow
  `inla.mode="classic"` to use proper line search.
- Avoid deprecated Matrix (\>=1.4-2) class coercion methods
- Work around for lack of full Matrix and ModelMatrix support for the
  `unique` method. Fixes
  [\#145](https://github.com/inlabru-org/inlabru/issues/145)

## inlabru 2.5.2

CRAN release: 2022-03-30

- More robust package checks
- More robust namespace and INLA availability checks
- Add package vignette with links to the website examples

## inlabru 2.5.1

- Revert to R language features compatible with R 4.0.5
- Use `strategy="gaussian"` during iterations.

## inlabru 2.5.0

CRAN release: 2022-03-21

### Features

- Add [`bru()`](https://inlabru-org.github.io/inlabru/reference/bru.md)
  timing information in `$bru_timings` and `$bru_iinla$timings`
- Add `SpatialPolygonsDataFrame` support to
  [`gg()`](https://inlabru-org.github.io/inlabru/reference/gg.md)
  methods
- Allow accessing `E` and `Ntrials` from `response_data` and `data`
  (further special arguments remain to be added)
- `deltaIC` improvements
- New transformation helper tools
  `bru_{forward/inverse}_transformation()`
- Experimental support for matrix and formula component inputs. E.g.
  with `~ name(~ -1 + a + b + a:b, model = "fixed")`, covariate fixed
  effect interaction specifications can be made. For formula input,
  [`MatrixModels::model.Matrix()`](https://rdrr.io/pkg/MatrixModels/man/model.Matrix.html)
  is called to construct matrix input that is then used as the A-matrix
  for fixed effects, one per column, added up to form the combined
  effect.
- Documentation and examples improvements

### Bug fixes

- Fix A-matrix construction for
  [`evaluate_model()`](https://inlabru-org.github.io/inlabru/reference/evaluate_model.md)
  for cases where the `inla_f` argument matters
- More efficient and robust mesh integration code
- Cleanup of environment handling for component lists

## inlabru 2.4.0

CRAN release: 2021-12-19

### Features

- Allow predictors to have different size than the input data. The
  `data` argument is now allowed to be a
  [`list()`](https://rdrr.io/r/base/list.html), and the new argument
  `response_data` allows separate specification of component inputs and
  response variables.
- Add `bru_mapper_collect` class for handling sequential collections of
  mappers, including collections where all but the first mapper is
  hidden from the [`INLA::f()`](https://rdrr.io/pkg/INLA/man/f.html)
  arguments `n` and `values`, as needed to support e.g. “bym2” models.
- Add `control.family` as a direct argument to
  [`like()`](https://inlabru-org.github.io/inlabru/reference/bru_obs.md).
  Gives a warning if a `control.family` argument is supplied to the the
  `options` argument of
  [`bru()`](https://inlabru-org.github.io/inlabru/reference/bru.md), but
  at least one likelihood has `control.family` information. (Issue
  [\#109](https://github.com/inlabru-org/inlabru/issues/109))

### Bugfixes

- Fix support for `SpatialPointsDataFrame` and `SpatialGridDataFrame`
  input to
  [`bru_fill_missing()`](https://inlabru-org.github.io/inlabru/reference/bru_fill_missing.md)
- Force explicit `model = "offset"` components instead of special
  options, to avoid interfering with the linearisation system (Issue
  [\#123](https://github.com/inlabru-org/inlabru/issues/123))
- Make the iterations more robust by resetting the internal INLA
  predictor states to initial value zero at each step

### Miscellaneous

- Rename the option `bru_method$stop_at_max_rel_deviation` to
  `bru_method$rel_tol`. Automatic conversion to the new name, but a
  warning is given.
- Add option `bru_method$max_step` to control the largest allowed line
  search scaling factor. See
  [`?bru_options`](https://inlabru-org.github.io/inlabru/reference/bru_options.md)
- New default option `bru_compress_cp` set to `TRUE` to compress the
  predictor expression for `family="cp"` to use a single element for the
  linear predictor sum.

## inlabru 2.3.1

CRAN release: 2021-03-22

- Documentation and dependency updates for CRAN compatibility
- See NEWS for version 2.3.0 for the major updates since version 2.1.13

## inlabru 2.3.0

CRAN release: 2021-03-16

### Breaking changes since version 2.1.13

- The model component argument `map` has been deprecated. Use `main` to
  specify the main component input,
  `~ elev(main = elevation, model = "rw2")`. Unlike the old `map`
  argument, `main` is the first one, so the shorter version
  `~ elev(elevation, model = "rw2")` also works.
- Intercept-like components should now have explicit inputs,
  e.g. `~ Intercept(1)` to avoid accidental confusion with other
  variables.
- The argument list for
  [`bru()`](https://inlabru-org.github.io/inlabru/reference/bru.md) has
  been simplified, so that all arguments except `components` and
  `options` must either be outputs from calls to
  [`like()`](https://inlabru-org.github.io/inlabru/reference/bru_obs.md),
  or arguments that can be sent to a single
  [`like()`](https://inlabru-org.github.io/inlabru/reference/bru_obs.md)
  call.
- The option setting system has been replaced with a more coherent
  system; see `?bru_options()` for details.
- The `samplers` and `domain` system for `lgcp` models is now stricter,
  and requires explicit `domain` definitions for all the point process
  dimensions. Alternatively, user-defined integration schemes can be
  supplied via the `ips` argument.

### New features since version 2.1.13

- The model component input arguments `main`, `group`, `replicate`, and
  `weights` can now take general R expressions using the data inputs.
  Special cases are detected: `SpatialPixels/GridDataFrame` objects are
  evaluated at spatial locations if the input data is a
  `SpatialPointsDataFrame` object. Functions are evaluated on the data
  object, e.g. `field(coordinates, model = spde)`
- The component arguments `mapper`, `group_mapper`, and
  `replicate_mapper` can be used for precise control of the mapping
  between inputs and latent variables. See
  [`?bru_mapper`](https://inlabru-org.github.io/inlabru/reference/bru_mapper.md)
  for more details. Mapper information is automatically extracted from
  [`INLA::inla.spde2.pcmatern()`](https://rdrr.io/pkg/INLA/man/inla.spde2.pcmatern.html)
  model objects.
- The R-INLA `weights` and `copy` features are now supported.
- The predictor expressions can access the data object directly via
  `.data.`
- If data from several rows can affect the same output row, the
  `allow_combine = TRUE` argument must be supplied to
  [`like()`](https://inlabru-org.github.io/inlabru/reference/bru_obs.md)
- The `include` and `exclude` arguments to
  [`like()`](https://inlabru-org.github.io/inlabru/reference/bru_obs.md),
  [`generate()`](https://inlabru-org.github.io/inlabru/reference/generate.md),
  and
  [`predict()`](https://rspatial.github.io/terra/reference/predict.html)
  can be used to specify which components are used for a given
  likelihood model or predictor expression. This can be used to prevent
  evaluation of components that are invalid for a likelihood or
  predictor.
- Predictor expressions can access the latent state of a model component
  directly, by adding the suffix `_latent` to the component name,
  e.g. `name_latent`. For
  [`like()`](https://inlabru-org.github.io/inlabru/reference/bru_obs.md),
  this requires `allow_latent = TRUE` to activate the needed
  linearisation code for this.
- Predictor expressions can evaluate component effects for arbitrary
  inputs by adding the suffix `_eval` to access special evaluator
  functions, e.g. `name_eval(1:10)`. This is useful for evaluating the
  1D effect of spatial covariates. See the NEWS item for version 2.2.8
  for further details.
- The internal system for predictor linearisation and iterated INLA
  inference has been rewritten to be faster and more robust
- See the NEWS entries for versions 2.1.14 to 2.2.8 for further details
  on new features and bug fixes

## inlabru 2.2.8

- Add `_eval` suffix feature for `generate.bru` and `predict.bru`, that
  provides a general evaluator function for each component, allowing
  evaluation of e.g. nonlinear effects of spatial covariates as a
  function of the covariate value instead of the by the spatial
  evaluator used in the component definition. For example, with
  `components = ~ covar(spatial_grid_df, model = "rw1")`, the prediction
  expression can have `~ covar_eval(covariate)`, where `covariate` is a
  data column in the prediction data object.

  For components with `group` and `replicate` features, these also need
  to be provided to the `_eval` function, with
  `..._eval(..., group = ..., replicate = ...)`

  This feature is built on top of the `_latent` suffix feature, that
  gives direct access to the latent state variables of a component, so
  in order to use `_eval` in the model predictor itself, you must use
  `like(..., allow_latent = TRUE)` in the model definition.

## inlabru 2.2.7

- Add support for `ngroup` and `nrep` in component definitions
- Updated `mexdolphin` and `mrsea` data sets, with consistent km units
  and improved mesh designs

## inlabru 2.2.6

- Add `predict(..., include)` discussion to distance sampling vignette,
  for handling non-spatial prediction in spatial models.
- Fix bugs in `gg.SpatialLines`

## inlabru 2.2.5

- Vignette corrections
- Documentation improvements
- Fix minor bug in `Spatial*` object handling and plotting

## inlabru 2.2.4

- Properly extract the joint latent conditional mode instead of the
  marginal latent conditional mode

## inlabru 2.2.2

- Fixed issue with
  [`predict()`](https://rspatial.github.io/terra/reference/predict.html)
  logic for converting output to `Spatial*DataFrame`
- Use `control.mode=list(restart=FALSE)` in the final inla run for
  nonlinear models, to avoid an unnecessary optimisation.
- Fix issues in `pixels()` and
  [`bru_fill_missing()`](https://inlabru-org.github.io/inlabru/reference/bru_fill_missing.md)
  for `Spatial*DataFrame` objects with `ncol=0` data frame parts.

## inlabru 2.2.1

- Fixed code regression bug for function input of covariates

## inlabru 2.2.0

- Support for the INLA “copy” feature, `comp2(input, copy = "comp1")`
- Allow component weights to be an unnamed parameter,
  `comp(input, weights, ...)`
- Direct access to the data objects in component inputs and predictor
  expressions, as `.data.`, allowing e.g. `covar(fun(.data.), ...)` for
  a complex covariate extractor method `fun()`
- Partial support for spherical manifold meshes
- Uses INLA integration strategy “eb” for initial nonlinear iterations,
  and a specified integration strategy only for the final iteration, so
  that the computations are faster, and uses the conditional latent mode
  as linearisation point.

## inlabru 2.1.15

- New options system
- New faster linearisation method
- New line search method to make the nonlinear inla iterations robust
- Method for updating old stored estimation objects
- System for supplying mappings between latent models and evaluated
  effects via `bru_mapper` objects
- Improved factor support; Either as “contrast with the 1st level”, via
  the special `"factor_contrast"` model, or all levels with model
  `"factor_full"`. Further options planned (e.g. a simpler options to
  fix the precision parameter). The estimated coefficients appear as
  random effects in the
  [`inla()`](https://rdrr.io/pkg/INLA/man/inla.html) output.
- Interface restructuring to support new features while keeping most
  backwards compatibility. Change `map=` to `main=` or unnamed first
  argument; Since `main` is the first parameter, it doesn’t need to be a
  named argument.
- Keep components with zero derivative in the linearisation
- PROJ6 support
- Add random seed option for posterior sampling
- Add package unit testing
- New backend code to make extended feature support easier
- New `int.args` option to control spatial integration resolution,
  thanks to Martin Jullum (`martinju`)

## inlabru 2.1.13

CRAN release: 2020-02-16

- Fix CRAN complaint regarding documentation

## inlabru 2.1.12

CRAN release: 2019-06-24

- Workaround an integration points error for old (ca pre-2018) INLA
  versions

## inlabru 2.1.11

- Add CITATION file

## inlabru 2.1.10

- Fix internal sampling bug due to INLA changes

## inlabru 2.1.9

CRAN release: 2018-07-24

- Remove unused `VignetteBuilder` entry from `DESCRIPTION`

## inlabru 2.1.8

- Update default options
- Prevent `int.polygon` from integrating outside the mesh domain, and
  generally more robust integration scheme construction.
- Fix [`bru()`](https://inlabru-org.github.io/inlabru/reference/bru.md)
  to
  [`like()`](https://inlabru-org.github.io/inlabru/reference/bru_obs.md)
  parameter logic. (Thanks to Peter Vesk for bug example)

## inlabru 2.1.7

- Added a `NEWS.md` file to track changes to the package.
- Added `inla` methods for
  [`predict()`](https://rspatial.github.io/terra/reference/predict.html)
  and
  [`generate()`](https://inlabru-org.github.io/inlabru/reference/generate.md)
  that convert `inla` output into `bru` objects before calling the `bru`
  prediction and posterior sample generator.
- Added protection for examples requiring optional packages
- Fix `sample.lgcp` output formatting, extended CRS support, and more
  efficient sampling algorithm
- Avoid dense matrices for effect mapping

## inlabru 2.1.4

- [`iinla()`](https://inlabru-org.github.io/inlabru/reference/iinla.md)
  tracks convergence of both fixed and random effects

## inlabru 2.1.3

CRAN release: 2018-02-11

- Added matrix geom
  [`gg.matrix()`](https://inlabru-org.github.io/inlabru/reference/gg.matrix.md)
- Fixed CRAN test issues
