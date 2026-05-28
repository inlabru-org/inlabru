# Information criteria for point processes

**\[experimental\]** Calculate information criteria for point process
models, using a fitted model and optionally new data. The interface is
still under development and may change without warning in future
versions.

## Usage

``` r
lgcp_IC(fit, data, predictor, domain, samplers, n.samples = 2000L)

delta_IC(data)
```

## Arguments

- fit:

  A fitted `bru` model object, from
  [`bru()`](https://inlabru-org.github.io/inlabru/reference/bru.md) or
  [`lgcp()`](https://inlabru-org.github.io/inlabru/reference/lgcp.md).

- data:

  A data frame containing the observed points, with columns for
  coordinates and any covariates used in the model. The data frame may
  also include a column named `.block` that indicates the block
  membership of each point, if blockwise criteria are to be calculated.

- predictor:

  An expression that evaluates to the linear predictor of the model,
  given the data. This should be an expression that can be evaluated in
  the context of the data frame and model components, just as in the
  `formula` argument of
  [`bru_obs()`](https://inlabru-org.github.io/inlabru/reference/bru_obs.md)
  and
  [`lgcp()`](https://inlabru-org.github.io/inlabru/reference/lgcp.md).

- domain, samplers:

  Arguments defining the integration domain(s) and sampling regions, as
  for
  [`bru_obs()`](https://inlabru-org.github.io/inlabru/reference/bru_obs.md)
  and
  [`lgcp()`](https://inlabru-org.github.io/inlabru/reference/lgcp.md).

- n.samples:

  Integer. The number of posterior samples to draw for estimating the
  criteria. Default is 2000.

## Value

A `tibble` with columns for the `Criterion` name, `Method`, log
pointwise predictive density (`lppd`), effective number of parameters
(`p_eff`), and the information criterion value (`IC`)

## Functions

- `delta_IC()`: Calculate IC differences based on joined output from
  several `lgcp_IC` calls;

      delta_IC(rbind(cbind(lgcp_IC(...), Model = "A"),
                     cbind(lgcp_IC(...), Model = "B")))
