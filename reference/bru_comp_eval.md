# Evaluate component values in predictor expressions

In predictor expressions, `name_eval(...)` can be used to evaluate the
effect of a component called "name".

## Usage

``` r
bru_comp_eval(
  main,
  group = NULL,
  replicate = NULL,
  weights = NULL,
  .state = NULL
)
```

## Arguments

- main, group, replicate, weights:

  Specification of where to evaluate a component. The four inputs are
  passed on to the joint `bru_mapper` for the component, as

      list(core = list(
             main = main,
             group = group,
             replicate = replicate),
           scale = weights)

  NOTE: If you have model component with the same name as a data
  variable you want to supply as input to `name_eval()`, you need to use
  `.data.[["myvar"]]` to access it. Otherwise, it will try to use the
  other component effect as input, which is ill-defined.

- .state:

  The internal component state. Normally supplied automatically by the
  internal methods for evaluating inlabru predictor expressions.

## Value

A vector of values for a component

## Examples

``` r
# \donttest{
if (bru_safe_inla() &&
  requireNamespace("sn", quietly = TRUE)) {
  set.seed(12345L)
  data <-
    data.frame(
      x = runif(5),
      idx = seq_len(5)
    )
  data$y <- data$x + rnorm(5, sd = 0.1)
  fit <- bru(
    y ~ 0 + x + field(idx, model = "iid", mapper = bm_index(5L)),
    family = "gaussian", data = data,
    control.family = list(
      hyper = list(prec = list(initial = 9, fixed = TRUE))
    ),
    options = list(control.inla = list(int.strategy = "eb"))
  )
  pred <- generate(
    fit,
    newdata = list(),
    formula = ~ field_eval(c(seq_len(5), 6, 7, 6)),
    n.samples = 1L
  )
}
# }
```
