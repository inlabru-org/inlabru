# Augment a bru model fit with fitted values

Adds posterior mean and credible interval columns to `data`. The user
must supply `pred_formula` because inlabru prediction expressions are
arbitrary R and cannot be recovered from the fit object alone. Also see
[`generics::augment()`](https://generics.r-lib.org/reference/augment.html).

## Usage

``` r
# S3 method for class 'bru'
augment(x, data, pred_formula, n_samples = 500L, seed = NULL, ...)
```

## Arguments

- x:

  A fitted `bru` object.

- data:

  A data frame of covariate values.

- pred_formula:

  A formula passed to
  [`inlabru::predict()`](https://rdrr.io/r/stats/predict.html).

- n_samples:

  Number of posterior samples (default 500).

- seed:

  Optional integer. If non-`NULL`, the session RNG is reset to this
  value before
  [`inlabru::predict()`](https://rdrr.io/r/stats/predict.html) is called
  so that two identical calls return identical posterior summaries.
  `NULL` leaves the session RNG state unchanged (predictions then vary
  run-to-run because
  [`inlabru::predict()`](https://rdrr.io/r/stats/predict.html) samples
  from the posterior).

- ...:

  Passed to
  [`inlabru::predict()`](https://rdrr.io/r/stats/predict.html).

## Value

`data` with additional columns `.fitted`, `.fitted_low`, `.fitted_high`,
`.fitted_sd`.
