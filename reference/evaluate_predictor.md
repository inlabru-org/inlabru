# Evaluate component effects or expressions

Evaluate component effects or expressions, based on a bru model and one
or several states of the latent variables and hyperparameters.

## Usage

``` r
evaluate_predictor(
  model,
  state,
  data,
  data_extra,
  effects,
  predictor,
  used = NULL,
  format = "auto",
  n_pred = NULL
)
```

## Arguments

- state:

  A list where each element is a list of named latent state information,
  as produced by
  [`evaluate_state()`](https://inlabru-org.github.io/inlabru/reference/evaluate_model.md)

- data:

  A `list`, `data.frame`, or `Spatial*DataFrame`, with coordinates and
  covariates needed to evaluate the model.

- data_extra:

  Additional data for the predictor evaluation. Variables with the same
  name as in `data` will be ignored, unless accessed via
  `.data_extra.[["name"]]` or `.data_extra.$name`, or via pronouns; see
  Details.

- effects:

  A list where each element is list of named evaluated effects, each
  computed by
  [`evaluate_effect_single_state.bru_comp_list()`](https://inlabru-org.github.io/inlabru/reference/evaluate_effect.md)

- predictor:

  Either a formula or
  [bru_pred_expr](https://inlabru-org.github.io/inlabru/reference/bru_pred_expr.md)
  expression

- used:

  A
  [`bru_used()`](https://inlabru-org.github.io/inlabru/reference/bru_used.md)
  object, or NULL (default)

- format:

  character; determines the storage format of the output. Available
  options:

  - `"auto"` If the first evaluated result is a vector or single-column
    matrix, the "matrix" format is used, otherwise "list".

  - `"matrix"` A matrix where each column contains the evaluated
    predictor expression for a state.

  - `"list"` A list where each column contains the evaluated predictor
    expression for a state.

  Default: "auto"

- n_pred:

  integer. If provided, scalar predictor results are expanded to vectors
  of length `n_pred`.

## Value

A list or matrix is returned, as specified by `format`

## Details

For each component, e.g. "name", the latent state values are available
as `name_latent`, and arbitrary evaluation can be done with
`name_eval(...)`, see
[`bru_comp_eval()`](https://inlabru-org.github.io/inlabru/reference/bru_comp_eval.md).

The evaluation supports several
[`rlang::as_data_pronoun()`](https://rlang.r-lib.org/reference/as_data_mask.html)
data masking pronouns, to access variables from different data sources,
and some of these also have corresponding full objects, with an appended
`.` in the name. The full objects can be passed as arguments to
functions.

- .effect/.effect.:

  refers to the `effects` vectors

- .latent/.latent.:

  refers to the latent state vectors

- .data/.data.:

  refers to the main `data` argument

- .data_extra/.data_extra.:

  refers to the `data_extra` argument

- .env:

  refers to the evaluation environment of the predictor
