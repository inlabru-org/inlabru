# Package index

## Model construction and estimation

Functions for constructing models and performing estimation.

- [`inlabru-package`](https://inlabru-org.github.io/inlabru/reference/inlabru-package.md)
  [`inlabru`](https://inlabru-org.github.io/inlabru/reference/inlabru-package.md)
  : inlabru

- [`bru()`](https://inlabru-org.github.io/inlabru/reference/bru.md)
  [`bru_rerun()`](https://inlabru-org.github.io/inlabru/reference/bru.md)
  [`summary(`*`<bru>`*`)`](https://inlabru-org.github.io/inlabru/reference/bru.md)
  [`print(`*`<summary_bru>`*`)`](https://inlabru-org.github.io/inlabru/reference/bru.md)
  [`print(`*`<bru>`*`)`](https://inlabru-org.github.io/inlabru/reference/bru.md)
  : Convenient model fitting using (iterated) INLA

- [`bru_obs()`](https://inlabru-org.github.io/inlabru/reference/bru_obs.md)
  [`like()`](https://inlabru-org.github.io/inlabru/reference/bru_obs.md)
  [`bru_obs_list()`](https://inlabru-org.github.io/inlabru/reference/bru_obs.md)
  [`c(`*`<bru_obs>`*`)`](https://inlabru-org.github.io/inlabru/reference/bru_obs.md)
  [`c(`*`<bru_obs_list>`*`)`](https://inlabru-org.github.io/inlabru/reference/bru_obs.md)
  [`` `[`( ``*`<bru_obs_list>`*`)`](https://inlabru-org.github.io/inlabru/reference/bru_obs.md)
  [`like_list()`](https://inlabru-org.github.io/inlabru/reference/bru_obs.md)
  [`bru_like_list()`](https://inlabru-org.github.io/inlabru/reference/bru_obs.md)
  :

  Observation model construction for usage with
  [`bru()`](https://inlabru-org.github.io/inlabru/reference/bru.md)

- [`bru_options()`](https://inlabru-org.github.io/inlabru/reference/bru_options.md)
  [`as.bru_options()`](https://inlabru-org.github.io/inlabru/reference/bru_options.md)
  [`bru_options_default()`](https://inlabru-org.github.io/inlabru/reference/bru_options.md)
  [`bru_options_check()`](https://inlabru-org.github.io/inlabru/reference/bru_options.md)
  [`bru_options_get()`](https://inlabru-org.github.io/inlabru/reference/bru_options.md)
  [`bru_options_set()`](https://inlabru-org.github.io/inlabru/reference/bru_options.md)
  [`bru_options_reset()`](https://inlabru-org.github.io/inlabru/reference/bru_options.md)
  [`bru_options_set_local()`](https://inlabru-org.github.io/inlabru/reference/bru_options.md)
  : Create or update an options objects

- [`bru_comp()`](https://inlabru-org.github.io/inlabru/reference/bru_comp.md)
  [`bru_component()`](https://inlabru-org.github.io/inlabru/reference/bru_comp.md)
  : Latent model component construction

- [`bru_comp_env_extra()`](https://inlabru-org.github.io/inlabru/reference/bru_comp_env_extra.md)
  [`` `bru_comp_env_extra<-`() ``](https://inlabru-org.github.io/inlabru/reference/bru_comp_env_extra.md)
  [`bru_comp_env()`](https://inlabru-org.github.io/inlabru/reference/bru_comp_env_extra.md)
  [`` `bru_comp_env<-`() ``](https://inlabru-org.github.io/inlabru/reference/bru_comp_env_extra.md)
  : Get/set component environment data

- [`bru_comp_eval()`](https://inlabru-org.github.io/inlabru/reference/bru_comp_eval.md)
  : Evaluate component values in predictor expressions

- [`bru_comp_list()`](https://inlabru-org.github.io/inlabru/reference/bru_comp_list.md)
  [`c(`*`<bru_comp_list>`*`)`](https://inlabru-org.github.io/inlabru/reference/bru_comp_list.md)
  [`c(`*`<bru_comp>`*`)`](https://inlabru-org.github.io/inlabru/reference/bru_comp_list.md)
  [`` `[`( ``*`<bru_comp_list>`*`)`](https://inlabru-org.github.io/inlabru/reference/bru_comp_list.md)
  : Methods for inlabru component lists

- [`new_bru_input()`](https://inlabru-org.github.io/inlabru/reference/bru_input.md)
  [`bru_input()`](https://inlabru-org.github.io/inlabru/reference/bru_input.md)
  : Store and evaluate component inputs

- [`lgcp()`](https://inlabru-org.github.io/inlabru/reference/lgcp.md) :
  Log Gaussian Cox process (LGCP) inference using INLA

- [`bru_response_size()`](https://inlabru-org.github.io/inlabru/reference/bru_response_size.md)
  : Response size queries

- [`bru_set_missing()`](https://inlabru-org.github.io/inlabru/reference/bru_set_missing.md)
  [`` `bru_set_missing<-`() ``](https://inlabru-org.github.io/inlabru/reference/bru_set_missing.md)
  : Set missing values in observation models

- [`bru_index()`](https://inlabru-org.github.io/inlabru/reference/bru_index.md)
  [`index_eval()`](https://inlabru-org.github.io/inlabru/reference/bru_index.md)
  **\[experimental\]** : Extract predictor or component index
  information

## Posterior evaluation

Functions for evaluating posterior properties.

- [`bru_names()`](https://inlabru-org.github.io/inlabru/reference/bru_names.md)
  : Extract standardised names from a bru or inla result object

- [`generate()`](https://inlabru-org.github.io/inlabru/reference/generate.md)
  : Generate samples from fitted bru models

- [`predict(`*`<bru>`*`)`](https://inlabru-org.github.io/inlabru/reference/predict.bru.md)
  : Prediction from fitted bru model

- [`spde.posterior()`](https://inlabru-org.github.io/inlabru/reference/spde.posterior.md)
  : Posteriors of SPDE hyper parameters and Matern correlation or
  covariance function.

- [`deltaIC()`](https://inlabru-org.github.io/inlabru/reference/deltaIC.md)
  :

  Summarise DIC from `lgcp` estimates.

- [`devel.cvmeasure()`](https://inlabru-org.github.io/inlabru/reference/devel.cvmeasure.md)
  : Variance and correlations measures for prediction components

## Optimization log information

Accessing the optimization text log, and plotting the optimization
convergence.

- [`bru_log()`](https://inlabru-org.github.io/inlabru/reference/bru_log.md)
  [`format(`*`<bru_log>`*`)`](https://inlabru-org.github.io/inlabru/reference/bru_log.md)
  [`print(`*`<bru_log>`*`)`](https://inlabru-org.github.io/inlabru/reference/bru_log.md)
  [`as.character(`*`<bru_log>`*`)`](https://inlabru-org.github.io/inlabru/reference/bru_log.md)
  [`` `[`( ``*`<bru_log>`*`)`](https://inlabru-org.github.io/inlabru/reference/bru_log.md)
  [`c(`*`<bru_log>`*`)`](https://inlabru-org.github.io/inlabru/reference/bru_log.md)
  [`length(`*`<bru_log>`*`)`](https://inlabru-org.github.io/inlabru/reference/bru_log.md)
  :

  Access methods for `bru_log` objects

- [`bru_log_bookmark()`](https://inlabru-org.github.io/inlabru/reference/bru_log_bookmark.md)
  [`bru_log_bookmarks()`](https://inlabru-org.github.io/inlabru/reference/bru_log_bookmark.md)
  :

  Methods for `bru_log` bookmarks

- [`bru_log_message()`](https://inlabru-org.github.io/inlabru/reference/bru_log_message.md)
  [`bru_log_abort()`](https://inlabru-org.github.io/inlabru/reference/bru_log_message.md)
  [`bru_log_warn()`](https://inlabru-org.github.io/inlabru/reference/bru_log_message.md)
  : Add a log message

- [`bru_log_new()`](https://inlabru-org.github.io/inlabru/reference/bru_log_new.md)
  :

  Create a `bru_log` object

- [`bru_log_offset()`](https://inlabru-org.github.io/inlabru/reference/bru_log_offset.md)
  [`bru_log_index()`](https://inlabru-org.github.io/inlabru/reference/bru_log_offset.md)
  :

  Position methods for `bru_log` objects

- [`bru_log_reset()`](https://inlabru-org.github.io/inlabru/reference/bru_log_reset.md)
  : Clear log contents

- [`bru_convergence_plot()`](https://inlabru-org.github.io/inlabru/reference/bru_convergence_plot.md)
  : Plot inlabru convergence diagnostics

- [`bru_timings()`](https://inlabru-org.github.io/inlabru/reference/bru_timings.md)
  : Extract timing information from fitted bru object

- [`bru_timings_plot()`](https://inlabru-org.github.io/inlabru/reference/bru_timings_plot.md)
  : Plot inlabru iteration timings

## Object/property extraction/queries

Functions for extracting or querying properties of objects.

- [`as_bru_comp()`](https://inlabru-org.github.io/inlabru/reference/as_bru_comp.md)
  [`as_bru_comp_list()`](https://inlabru-org.github.io/inlabru/reference/as_bru_comp.md)
  :

  Conversion methods for `bru_comp` and `bru_comp_list` objects

- [`as_bru_mapper()`](https://inlabru-org.github.io/inlabru/reference/as_bru_mapper.md)
  : Methods for mapper extraction

- [`as_bru_obs()`](https://inlabru-org.github.io/inlabru/reference/as_bru_obs.md)
  [`as_bru_obs_list()`](https://inlabru-org.github.io/inlabru/reference/as_bru_obs.md)
  :

  Conversion methods for `bru_obs` and `bru_obs_list` objects

- [`bru_info()`](https://inlabru-org.github.io/inlabru/reference/bru_info.md)
  [`as_bru_info()`](https://inlabru-org.github.io/inlabru/reference/bru_info.md)
  [`summary(`*`<bru_info>`*`)`](https://inlabru-org.github.io/inlabru/reference/bru_info.md)
  [`print(`*`<summary_bru_info>`*`)`](https://inlabru-org.github.io/inlabru/reference/bru_info.md)
  [`print(`*`<bru_info>`*`)`](https://inlabru-org.github.io/inlabru/reference/bru_info.md)
  : Methods for bru_info objects

- [`bru_is_additive()`](https://inlabru-org.github.io/inlabru/reference/bru_is_additive.md)
  : Check for predictor expression additivity

- [`bru_is_linear()`](https://inlabru-org.github.io/inlabru/reference/bru_is_linear.md)
  : Check for predictor linearity

- [`bru_is_rowwise()`](https://inlabru-org.github.io/inlabru/reference/bru_is_rowwise.md)
  : Check for predictor rowwise evaluability

## Mapper methods

Functions for `bru_mapper` handling.

### Low level and wrapper methods

- [`bru_get_mapper()`](https://inlabru-org.github.io/inlabru/reference/bru_get_mapper.md)
  [`bru_get_mapper_safely()`](https://inlabru-org.github.io/inlabru/reference/bru_get_mapper.md)
  : Extract mapper information from INLA model component objects

- [`bru_mapper()`](https://inlabru-org.github.io/inlabru/reference/bru_mapper.md)
  [`bru_mapper_define()`](https://inlabru-org.github.io/inlabru/reference/bru_mapper.md)
  :

  Constructors for `bru_mapper` objects

### Constructor methods

- [`bm_aggregate()`](https://inlabru-org.github.io/inlabru/reference/bm_aggregate.md)
  [`bru_mapper_aggregate()`](https://inlabru-org.github.io/inlabru/reference/bm_aggregate.md)
  : Mapper for aggregation

- [`bm_collect()`](https://inlabru-org.github.io/inlabru/reference/bm_collect.md)
  [`bru_mapper_collect()`](https://inlabru-org.github.io/inlabru/reference/bm_collect.md)
  [`` `[`( ``*`<bm_collect>`*`)`](https://inlabru-org.github.io/inlabru/reference/bm_collect.md)
  [`` `[`( ``*`<bru_mapper_collect>`*`)`](https://inlabru-org.github.io/inlabru/reference/bm_collect.md)
  : Mapper for concatenated variables

- [`bm_const()`](https://inlabru-org.github.io/inlabru/reference/bm_const.md)
  [`bru_mapper_const()`](https://inlabru-org.github.io/inlabru/reference/bm_const.md)
  : Constant mapper

- [`bm_factor()`](https://inlabru-org.github.io/inlabru/reference/bm_factor.md)
  [`bru_mapper_factor()`](https://inlabru-org.github.io/inlabru/reference/bm_factor.md)
  : Mapper for factor variables

- [`bru_mapper(`*`<fm_mesh_1d>`*`)`](https://inlabru-org.github.io/inlabru/reference/bm_fm_mesh_1d.md)
  :

  Mapper for `fm_mesh_1d`

- [`bm_fmesher()`](https://inlabru-org.github.io/inlabru/reference/bm_fmesher.md)
  [`bru_mapper_fmesher()`](https://inlabru-org.github.io/inlabru/reference/bm_fmesher.md)
  [`bru_mapper(`*`<fm_mesh_2d>`*`)`](https://inlabru-org.github.io/inlabru/reference/bm_fmesher.md)
  :

  Mapper for general `fmesher` function space objects

- [`bm_harmonics()`](https://inlabru-org.github.io/inlabru/reference/bm_harmonics.md)
  [`bru_mapper_harmonics()`](https://inlabru-org.github.io/inlabru/reference/bm_harmonics.md)
  : Mapper for cos/sin functions

- [`bm_index()`](https://inlabru-org.github.io/inlabru/reference/bm_index.md)
  [`bru_mapper_index()`](https://inlabru-org.github.io/inlabru/reference/bm_index.md)
  : Mapper for indexed variables

- [`bm_linear()`](https://inlabru-org.github.io/inlabru/reference/bm_linear.md)
  [`bru_mapper_linear()`](https://inlabru-org.github.io/inlabru/reference/bm_linear.md)
  : Mapper for a linear effect

- [`as_bm_list()`](https://inlabru-org.github.io/inlabru/reference/bm_list.md)
  [`c(`*`<bru_mapper>`*`)`](https://inlabru-org.github.io/inlabru/reference/bm_list.md)
  [`c(`*`<bm_list>`*`)`](https://inlabru-org.github.io/inlabru/reference/bm_list.md)
  [`` `[`( ``*`<bm_list>`*`)`](https://inlabru-org.github.io/inlabru/reference/bm_list.md)
  : Methods for mapper lists

- [`bm_logitaverage()`](https://inlabru-org.github.io/inlabru/reference/bm_logitaverage.md)
  **\[experimental\]** : Mapper for logit-sum-inverse-logit aggregation

- [`bm_logsumexp()`](https://inlabru-org.github.io/inlabru/reference/bm_logsumexp.md)
  [`bru_mapper_logsumexp()`](https://inlabru-org.github.io/inlabru/reference/bm_logsumexp.md)
  : Mapper for log-sum-exp aggregation

- [`bm_marginal()`](https://inlabru-org.github.io/inlabru/reference/bm_marginal.md)
  [`bru_mapper_marginal()`](https://inlabru-org.github.io/inlabru/reference/bm_marginal.md)
  : Mapper for marginal distribution transformation

- [`bm_matrix()`](https://inlabru-org.github.io/inlabru/reference/bm_matrix.md)
  [`bru_mapper_matrix()`](https://inlabru-org.github.io/inlabru/reference/bm_matrix.md)
  : Mapper for matrix multiplication

- [`bm_multi()`](https://inlabru-org.github.io/inlabru/reference/bm_multi.md)
  [`bru_mapper_multi()`](https://inlabru-org.github.io/inlabru/reference/bm_multi.md)
  [`` `[`( ``*`<bm_multi>`*`)`](https://inlabru-org.github.io/inlabru/reference/bm_multi.md)
  [`` `[`( ``*`<bru_mapper_multi>`*`)`](https://inlabru-org.github.io/inlabru/reference/bm_multi.md)
  : Mapper for tensor product domains

- [`bm_pipe()`](https://inlabru-org.github.io/inlabru/reference/bm_pipe.md)
  [`bru_mapper_pipe()`](https://inlabru-org.github.io/inlabru/reference/bm_pipe.md)
  : Mapper for linking several mappers in sequence

- [`bm_reparam()`](https://inlabru-org.github.io/inlabru/reference/bm_reparam.md)
  : Mapper for reparameterising mapper states

- [`bm_repeat()`](https://inlabru-org.github.io/inlabru/reference/bm_repeat.md)
  [`bru_mapper_repeat()`](https://inlabru-org.github.io/inlabru/reference/bm_repeat.md)
  : Mapper for repeating a mapper

- [`bm_scale()`](https://inlabru-org.github.io/inlabru/reference/bm_scale.md)
  [`bru_mapper_scale()`](https://inlabru-org.github.io/inlabru/reference/bm_scale.md)
  : Mapper for element-wise scaling

- [`bm_shift()`](https://inlabru-org.github.io/inlabru/reference/bm_shift.md)
  [`bru_mapper_shift()`](https://inlabru-org.github.io/inlabru/reference/bm_shift.md)
  : Mapper for element-wise shifting

- [`bm_sum()`](https://inlabru-org.github.io/inlabru/reference/bm_sum.md)
  [`bru_mapper_sum()`](https://inlabru-org.github.io/inlabru/reference/bm_sum.md)
  [`` `[`( ``*`<bm_sum>`*`)`](https://inlabru-org.github.io/inlabru/reference/bm_sum.md)
  [`` `[`( ``*`<bru_mapper_sum>`*`)`](https://inlabru-org.github.io/inlabru/reference/bm_sum.md)
  : Mapper for adding multiple mappers

- [`bm_taylor()`](https://inlabru-org.github.io/inlabru/reference/bm_taylor.md)
  [`bru_mapper_taylor()`](https://inlabru-org.github.io/inlabru/reference/bm_taylor.md)
  : Mapper for linear Taylor approximations

- [`ibm_input_set()`](https://inlabru-org.github.io/inlabru/reference/ibm_input.md)
  [`ibm_input_new()`](https://inlabru-org.github.io/inlabru/reference/ibm_input.md)
  [`ibm_input_available()`](https://inlabru-org.github.io/inlabru/reference/ibm_input.md)
  [`ibm_input_get()`](https://inlabru-org.github.io/inlabru/reference/ibm_input.md)
  [`bm_autodetect()`](https://inlabru-org.github.io/inlabru/reference/ibm_input.md)
  :

  Interface between `bru_input` and `bru_mapper`

### Evaluator methods

- [`bru_mapper_generics`](https://inlabru-org.github.io/inlabru/reference/bru_mapper_generics.md)
  : Generic methods for bru_mapper objects

- [`ibm_linear(`*`<bru_model>`*`)`](https://inlabru-org.github.io/inlabru/reference/bru_model_mapper_methods.md)
  [`ibm_linear(`*`<bru_comp_list>`*`)`](https://inlabru-org.github.io/inlabru/reference/bru_model_mapper_methods.md)
  [`ibm_simplify(`*`<bru_model>`*`)`](https://inlabru-org.github.io/inlabru/reference/bru_model_mapper_methods.md)
  [`ibm_simplify(`*`<bru_comp>`*`)`](https://inlabru-org.github.io/inlabru/reference/bru_model_mapper_methods.md)
  [`ibm_simplify(`*`<bru_comp_list>`*`)`](https://inlabru-org.github.io/inlabru/reference/bru_model_mapper_methods.md)
  [`ibm_linear(`*`<bm_list>`*`)`](https://inlabru-org.github.io/inlabru/reference/bru_model_mapper_methods.md)
  [`ibm_simplify(`*`<bm_list>`*`)`](https://inlabru-org.github.io/inlabru/reference/bru_model_mapper_methods.md)
  : Mapper methods for model objects

- [`ibm_eval()`](https://inlabru-org.github.io/inlabru/reference/ibm_eval.md)
  : Evaluate a mapping

- [`ibm_eval2()`](https://inlabru-org.github.io/inlabru/reference/ibm_eval2.md)
  : Evaluate a mapper and its Jacobian

- [`ibm_inla_subset()`](https://inlabru-org.github.io/inlabru/reference/ibm_inla_subset.md)
  : Find index subset of INLA visible states

- [`ibm_input_set()`](https://inlabru-org.github.io/inlabru/reference/ibm_input.md)
  [`ibm_input_new()`](https://inlabru-org.github.io/inlabru/reference/ibm_input.md)
  [`ibm_input_available()`](https://inlabru-org.github.io/inlabru/reference/ibm_input.md)
  [`ibm_input_get()`](https://inlabru-org.github.io/inlabru/reference/ibm_input.md)
  [`bm_autodetect()`](https://inlabru-org.github.io/inlabru/reference/ibm_input.md)
  :

  Interface between `bru_input` and `bru_mapper`

- [`ibm_invalid_output()`](https://inlabru-org.github.io/inlabru/reference/ibm_invalid_output.md)
  : Detect invalid input to a mapper

- [`ibm_is_linear()`](https://inlabru-org.github.io/inlabru/reference/ibm_is_linear.md)
  : Check if a mapper is linear/affine

- [`ibm_is_rowwise()`](https://inlabru-org.github.io/inlabru/reference/ibm_is_rowwise.md)
  : Check if a mapper is rowwise

- [`ibm_jacobian()`](https://inlabru-org.github.io/inlabru/reference/ibm_jacobian.md)
  : Jacobian of a mapper

- [`ibm_linear()`](https://inlabru-org.github.io/inlabru/reference/ibm_linear.md)
  : Compute a mapper linearisation

- [`ibm_n()`](https://inlabru-org.github.io/inlabru/reference/ibm_n.md)
  : Size of the latent vector of a mapping

- [`ibm_n_output()`](https://inlabru-org.github.io/inlabru/reference/ibm_n_output.md)
  : Output size of a mapping

- [`ibm_names()`](https://inlabru-org.github.io/inlabru/reference/ibm_names.md)
  [`` `ibm_names<-`() ``](https://inlabru-org.github.io/inlabru/reference/ibm_names.md)
  : Names of submapper

- [`ibm_simplify()`](https://inlabru-org.github.io/inlabru/reference/ibm_simplify.md)
  : Simplify a mapper

- [`ibm_values()`](https://inlabru-org.github.io/inlabru/reference/ibm_values.md)
  : Value vector for a mapping

### Conversion and helper methods

- [`as_bru_mapper()`](https://inlabru-org.github.io/inlabru/reference/as_bru_mapper.md)
  : Methods for mapper extraction
- [`as_bm_list()`](https://inlabru-org.github.io/inlabru/reference/bm_list.md)
  [`c(`*`<bru_mapper>`*`)`](https://inlabru-org.github.io/inlabru/reference/bm_list.md)
  [`c(`*`<bm_list>`*`)`](https://inlabru-org.github.io/inlabru/reference/bm_list.md)
  [`` `[`( ``*`<bm_list>`*`)`](https://inlabru-org.github.io/inlabru/reference/bm_list.md)
  : Methods for mapper lists
- [`bru_forward_transformation()`](https://inlabru-org.github.io/inlabru/reference/bru_transformation.md)
  [`bru_inverse_transformation()`](https://inlabru-org.github.io/inlabru/reference/bru_transformation.md)
  : Transformation tools

## Miscellaneous

Miscellaneous functions.

- [`eval_spatial()`](https://inlabru-org.github.io/inlabru/reference/eval_spatial.md)
  : Evaluate spatial covariates
- [`bru_fill_missing()`](https://inlabru-org.github.io/inlabru/reference/bru_fill_missing.md)
  : Fill in missing values in Spatial grids
- [`point2count()`](https://inlabru-org.github.io/inlabru/reference/point2count.md)
  : Convert a plot sample of points into one of counts.
- [`sample.lgcp()`](https://inlabru-org.github.io/inlabru/reference/sample.lgcp.md)
  : Sample from an inhomogeneous Poisson process

## Plot helpers

Functions for helping build figures.

- [`gg(`*`<data.frame>`*`)`](https://inlabru-org.github.io/inlabru/reference/gg.bru_prediction.md)
  [`gg(`*`<bru_prediction>`*`)`](https://inlabru-org.github.io/inlabru/reference/gg.bru_prediction.md)
  [`gg(`*`<prediction>`*`)`](https://inlabru-org.github.io/inlabru/reference/gg.bru_prediction.md)
  [`plot(`*`<bru_prediction>`*`)`](https://inlabru-org.github.io/inlabru/reference/gg.bru_prediction.md)
  [`plot(`*`<prediction>`*`)`](https://inlabru-org.github.io/inlabru/reference/gg.bru_prediction.md)
  : Geom for predictions
- [`plot(`*`<bru>`*`)`](https://inlabru-org.github.io/inlabru/reference/plot.bru.md)
  [`plotmarginal.inla()`](https://inlabru-org.github.io/inlabru/reference/plot.bru.md)
  : Plot method for posterior marginals estimated by bru
- [`plotsample()`](https://inlabru-org.github.io/inlabru/reference/plotsample.md)
  : Create a plot sample.
- [`glplot()`](https://inlabru-org.github.io/inlabru/reference/glplot.md)
  [`globe()`](https://inlabru-org.github.io/inlabru/reference/glplot.md)
  : Render objects using RGL
- [`gg(`*`<RasterLayer>`*`)`](https://inlabru-org.github.io/inlabru/reference/gg.RasterLayer.md)
  : Geom for RasterLayer objects
- [`gg()`](https://inlabru-org.github.io/inlabru/reference/gg.md) :
  ggplot2 geomes for inlabru related objects
- [`gg(`*`<SpatRaster>`*`)`](https://inlabru-org.github.io/inlabru/reference/gg.SpatRaster.md)
  : Geom wrapper for SpatRaster objects
- [`gg(`*`<SpatialPoints>`*`)`](https://inlabru-org.github.io/inlabru/reference/gg.Spatial.md)
  [`gg(`*`<SpatialLines>`*`)`](https://inlabru-org.github.io/inlabru/reference/gg.Spatial.md)
  [`gg(`*`<SpatialPolygons>`*`)`](https://inlabru-org.github.io/inlabru/reference/gg.Spatial.md)
  [`gg(`*`<SpatialGridDataFrame>`*`)`](https://inlabru-org.github.io/inlabru/reference/gg.Spatial.md)
  [`gg(`*`<SpatialPixelsDataFrame>`*`)`](https://inlabru-org.github.io/inlabru/reference/gg.Spatial.md)
  [`gg(`*`<SpatialPixels>`*`)`](https://inlabru-org.github.io/inlabru/reference/gg.Spatial.md)
  : Geoms for sp Spatial objects
- [`gg(`*`<fm_mesh_1d>`*`)`](https://inlabru-org.github.io/inlabru/reference/gg.fm_mesh_1d.md)
  : Geom for fm_mesh_1d objects
- [`gg(`*`<fm_mesh_2d>`*`)`](https://inlabru-org.github.io/inlabru/reference/gg.fm_mesh_2d.md)
  : Geom for fm_mesh_2d objects
- [`gg(`*`<matrix>`*`)`](https://inlabru-org.github.io/inlabru/reference/gg.matrix.md)
  : Geom for matrix
- [`gg(`*`<sf>`*`)`](https://inlabru-org.github.io/inlabru/reference/gg.sf.md)
  : Geom helper for sf objects
- [`bincount()`](https://inlabru-org.github.io/inlabru/reference/bincount.md)
  : 1D LGCP bin count simulation and comparison with data

## Printing methods

Functions for printing.

- [`format(`*`<bru_mapper>`*`)`](https://inlabru-org.github.io/inlabru/reference/bm_summary.md)
  [`format(`*`<bm_list>`*`)`](https://inlabru-org.github.io/inlabru/reference/bm_summary.md)
  [`summary(`*`<bru_mapper>`*`)`](https://inlabru-org.github.io/inlabru/reference/bm_summary.md)
  [`format(`*`<bm_multi>`*`)`](https://inlabru-org.github.io/inlabru/reference/bm_summary.md)
  [`format(`*`<bm_pipe>`*`)`](https://inlabru-org.github.io/inlabru/reference/bm_summary.md)
  [`format(`*`<bm_collect>`*`)`](https://inlabru-org.github.io/inlabru/reference/bm_summary.md)
  [`format(`*`<bm_sum>`*`)`](https://inlabru-org.github.io/inlabru/reference/bm_summary.md)
  [`format(`*`<bm_repeat>`*`)`](https://inlabru-org.github.io/inlabru/reference/bm_summary.md)
  [`format(`*`<bm_reparam>`*`)`](https://inlabru-org.github.io/inlabru/reference/bm_summary.md)
  [`print(`*`<summary_bru_mapper>`*`)`](https://inlabru-org.github.io/inlabru/reference/bm_summary.md)
  [`print(`*`<bru_mapper>`*`)`](https://inlabru-org.github.io/inlabru/reference/bm_summary.md)
  [`print(`*`<bm_list>`*`)`](https://inlabru-org.github.io/inlabru/reference/bm_summary.md)
  : mapper object summaries

- [`bru()`](https://inlabru-org.github.io/inlabru/reference/bru.md)
  [`bru_rerun()`](https://inlabru-org.github.io/inlabru/reference/bru.md)
  [`summary(`*`<bru>`*`)`](https://inlabru-org.github.io/inlabru/reference/bru.md)
  [`print(`*`<summary_bru>`*`)`](https://inlabru-org.github.io/inlabru/reference/bru.md)
  [`print(`*`<bru>`*`)`](https://inlabru-org.github.io/inlabru/reference/bru.md)
  : Convenient model fitting using (iterated) INLA

- [`bru_info()`](https://inlabru-org.github.io/inlabru/reference/bru_info.md)
  [`as_bru_info()`](https://inlabru-org.github.io/inlabru/reference/bru_info.md)
  [`summary(`*`<bru_info>`*`)`](https://inlabru-org.github.io/inlabru/reference/bru_info.md)
  [`print(`*`<summary_bru_info>`*`)`](https://inlabru-org.github.io/inlabru/reference/bru_info.md)
  [`print(`*`<bru_info>`*`)`](https://inlabru-org.github.io/inlabru/reference/bru_info.md)
  : Methods for bru_info objects

- [`bru_log()`](https://inlabru-org.github.io/inlabru/reference/bru_log.md)
  [`format(`*`<bru_log>`*`)`](https://inlabru-org.github.io/inlabru/reference/bru_log.md)
  [`print(`*`<bru_log>`*`)`](https://inlabru-org.github.io/inlabru/reference/bru_log.md)
  [`as.character(`*`<bru_log>`*`)`](https://inlabru-org.github.io/inlabru/reference/bru_log.md)
  [`` `[`( ``*`<bru_log>`*`)`](https://inlabru-org.github.io/inlabru/reference/bru_log.md)
  [`c(`*`<bru_log>`*`)`](https://inlabru-org.github.io/inlabru/reference/bru_log.md)
  [`length(`*`<bru_log>`*`)`](https://inlabru-org.github.io/inlabru/reference/bru_log.md)
  :

  Access methods for `bru_log` objects

- [`summary(`*`<bru_obs>`*`)`](https://inlabru-org.github.io/inlabru/reference/bru_obs_print.md)
  [`summary(`*`<bru_obs_list>`*`)`](https://inlabru-org.github.io/inlabru/reference/bru_obs_print.md)
  [`print(`*`<summary_bru_obs>`*`)`](https://inlabru-org.github.io/inlabru/reference/bru_obs_print.md)
  [`print(`*`<summary_bru_obs_list>`*`)`](https://inlabru-org.github.io/inlabru/reference/bru_obs_print.md)
  [`print(`*`<bru_obs>`*`)`](https://inlabru-org.github.io/inlabru/reference/bru_obs_print.md)
  [`print(`*`<bru_obs_list>`*`)`](https://inlabru-org.github.io/inlabru/reference/bru_obs_print.md)
  : Summary and print methods for observation models

- [`format(`*`<bru_input>`*`)`](https://inlabru-org.github.io/inlabru/reference/summary.bru_input.md)
  [`summary(`*`<bru_input>`*`)`](https://inlabru-org.github.io/inlabru/reference/summary.bru_input.md)
  [`print(`*`<bru_input>`*`)`](https://inlabru-org.github.io/inlabru/reference/summary.bru_input.md)
  : Summarise component inputs

- [`summary(`*`<bru_options>`*`)`](https://inlabru-org.github.io/inlabru/reference/summary.bru_options.md)
  [`print(`*`<summary_bru_options>`*`)`](https://inlabru-org.github.io/inlabru/reference/summary.bru_options.md)
  : Print inlabru options

## Data objects

Data objects

- [`Poisson1_1D`](https://inlabru-org.github.io/inlabru/reference/Poisson1_1D.md)
  : 1-Dimensional Homogeneous Poisson example.
- [`Poisson2_1D`](https://inlabru-org.github.io/inlabru/reference/Poisson2_1D.md)
  : 1-Dimensional NonHomogeneous Poisson example.
- [`Poisson3_1D`](https://inlabru-org.github.io/inlabru/reference/Poisson3_1D.md)
  : 1-Dimensional NonHomogeneous Poisson example.
- [`gorillas_sf`](https://inlabru-org.github.io/inlabru/reference/gorillas_sf.md)
  [`gorillas_sf_gcov()`](https://inlabru-org.github.io/inlabru/reference/gorillas_sf.md)
  [`gorillas_sp()`](https://inlabru-org.github.io/inlabru/reference/gorillas_sf.md)
  : Gorilla nesting sites in sf format
- [`mexdolphin_sf`](https://inlabru-org.github.io/inlabru/reference/mexdolphin_sf.md)
  [`mexdolphin_sp()`](https://inlabru-org.github.io/inlabru/reference/mexdolphin_sf.md)
  : Pan-tropical spotted dolphins in the Gulf of Mexico
- [`mrsea`](https://inlabru-org.github.io/inlabru/reference/mrsea.md) :
  Marine renewables strategic environmental assessment
- [`robins_subset`](https://inlabru-org.github.io/inlabru/reference/robins_subset.md)
  : robins_subset
- [`shrimp`](https://inlabru-org.github.io/inlabru/reference/shrimp.md)
  : Blue and red shrimp in the Western Mediterranean Sea
- [`toygroups`](https://inlabru-org.github.io/inlabru/reference/toygroups.md)
  : Simulated 1D animal group locations and group sizes
- [`toypoints`](https://inlabru-org.github.io/inlabru/reference/toypoints.md)
  : Simulated 2D point process data

## Internal helper functions

Helper functions for more low level operations.

- [`evaluate_effect_single_state()`](https://inlabru-org.github.io/inlabru/reference/evaluate_effect.md)
  : Evaluate a component effect

- [`evaluate_model()`](https://inlabru-org.github.io/inlabru/reference/evaluate_model.md)
  [`evaluate_state()`](https://inlabru-org.github.io/inlabru/reference/evaluate_model.md)
  : Evaluate or sample from a posterior result given a model and
  locations

- [`evaluate_predictor()`](https://inlabru-org.github.io/inlabru/reference/evaluate_predictor.md)
  : Evaluate component effects or expressions

- [`input_eval()`](https://inlabru-org.github.io/inlabru/reference/inlabru-deprecated.md)
  [`evaluate_inputs()`](https://inlabru-org.github.io/inlabru/reference/inlabru-deprecated.md)
  [`gmap()`](https://inlabru-org.github.io/inlabru/reference/inlabru-deprecated.md)
  [`gm()`](https://inlabru-org.github.io/inlabru/reference/inlabru-deprecated.md)
  [`row_kron()`](https://inlabru-org.github.io/inlabru/reference/inlabru-deprecated.md)
  : Deprecated functions in inlabru

- [`expand_labels()`](https://inlabru-org.github.io/inlabru/reference/expand_labels.md)
  : Expand labels

- [`bru_inla.stack.mjoin()`](https://inlabru-org.github.io/inlabru/reference/bru_inla.stack.mjoin.md)
  : Join stacks intended to be run with different likelihoods

- [`spatial.to.ppp()`](https://inlabru-org.github.io/inlabru/reference/spatial.to.ppp.md)
  : Convert SpatialPoints and boundary polygon to spatstat ppp object

- [`bru_make_stack()`](https://inlabru-org.github.io/inlabru/reference/bru_make_stack.md)
  : Build an inla data stack from linearisation information

- [`bru_summarise()`](https://inlabru-org.github.io/inlabru/reference/bru_summarise.md)
  : Summarise and annotate data

- [`bru_standardise_names()`](https://inlabru-org.github.io/inlabru/reference/bru_standardise_names.md)
  : Standardise inla hyperparameter names

- [`bru_safe_inla()`](https://inlabru-org.github.io/inlabru/reference/bru_safe_inla.md)
  : Load INLA safely for examples and tests

- [`bru_safe_sp()`](https://inlabru-org.github.io/inlabru/reference/bru_safe_sp.md)
  :

  Check for potential `sp` version compatibility issues

- [`bru_call_options()`](https://inlabru-org.github.io/inlabru/reference/bru_call_options.md)
  : Additional bru options

- [`bru_compute_linearisation()`](https://inlabru-org.github.io/inlabru/reference/bru_compute_linearisation.md)
  : Compute inlabru model linearisation information

- [`bru_is_additive()`](https://inlabru-org.github.io/inlabru/reference/bru_is_additive.md)
  : Check for predictor expression additivity

## Deprecated methods

Deprecated methods

- [`multiplot()`](https://inlabru-org.github.io/inlabru/reference/multiplot.md)
  **\[deprecated\]** : Multiple ggplots on a page.
