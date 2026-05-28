# Devel: Model evaluation flowchart

(Vignette under construction!)

## Mapping from component inputs and latent states to component effects

## Linearising a mapping

## Component input evaluation

For each `<label>` of `main`, `group`, `replicate`, and `weights`, the
given expression `expr` is evaluated in the data context, producing the
`input` to the component `mapper`. For spatial covariate inputs, the
corresponding `<label>_layer` expression is also evaluated.

Red nodes indicate deprecated behaviour retained for backwards
compatibility.

## Intergration point construction

Flow diagram for new integration scheme construction, implemented as
`fm_int(domain, samplers)` methods.
