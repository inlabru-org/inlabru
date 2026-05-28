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

      list(mapper = list(
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
if (bru_safe_inla() &&
  require("sf", quietly = TRUE) &&
  requireNamespace("sn", quietly = TRUE)) {
  mesh <- fmesher::fm_mesh_2d_inla(
    cbind(0, 0),
    offset = 2,
    max.edge = 2.5
  )
  spde <- INLA::inla.spde2.pcmatern(
    mesh,
    prior.range = c(1, NA),
    prior.sigma = c(0.2, NA)
  )
  set.seed(12345L)
  data <- sf::st_as_sf(
    data.frame(
      x = runif(50),
      y = runif(50),
      z = rnorm(50)
    ),
    coords = c("x", "y")
  )
  fit <- bru(
    z ~ -1 + field(geometry, model = spde),
    family = "gaussian", data = data,
    options = list(control.inla = list(int.strategy = "eb"))
  )
  pred <- generate(
    fit,
    newdata = data.frame(A = 0.5, B = 0.5),
    formula = ~ field_eval(cbind(A, B)),
    n.samples = 1L
  )
}
#> Linking to GEOS 3.12.1, GDAL 3.8.4, PROJ 9.4.0; sf_use_s2() is TRUE
```
