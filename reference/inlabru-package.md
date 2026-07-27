# inlabru

Convenient model fitting using (iterated) INLA.

## Details

`inlabru` facilitates Bayesian spatial modelling using integrated nested
Laplace approximations. It is heavily based on R-inla
(<https://www.r-inla.org>) but adds additional modelling abilities and
simplified syntax for (in particular) spatial models. Tutorials and more
information can be found at <https://inlabru-org.github.io/inlabru/> and
<http://www.inlabru.org/>. The iterative method used for non-linear
predictors is documented in the `method` vignette.

The main function for inference using inlabru is
[`bru()`](https://inlabru-org.github.io/inlabru/reference/bru.md). The
general model specification details is documented in
[`bru_comp()`](https://inlabru-org.github.io/inlabru/reference/bru_comp.md)
and
[`bru_obs()`](https://inlabru-org.github.io/inlabru/reference/bru_obs.md).
Posterior quantities beyond the basic summaries can be calculated with a
[`predict()`](https://rdrr.io/r/stats/predict.html) method, documented
in
[`predict.bru()`](https://inlabru-org.github.io/inlabru/reference/predict.md).
For point process inference
[`lgcp()`](https://inlabru-org.github.io/inlabru/reference/lgcp.md) can
be used as a shortcut to `bru(..., bru_obs(model="cp", ...))`.

The package comes with multiple real world data sets, namely
[gorillas](https://inlabru-org.github.io/inlabru/reference/gorillas.md),
[gorillas_sf](https://inlabru-org.github.io/inlabru/reference/gorillas_sf.md),
[mexdolphin_sf](https://inlabru-org.github.io/inlabru/reference/mexdolphin_sf.md).
Plotting these data sets is straight forward using inlabru's extensions
to `ggplot2`, e.g. the
[`gg()`](https://inlabru-org.github.io/inlabru/reference/gg.md)
function. For educational purposes some simulated data sets are
available as well, e.g.
[Poisson1_1D](https://inlabru-org.github.io/inlabru/reference/Poisson1_1D.md),
[Poisson2_1D](https://inlabru-org.github.io/inlabru/reference/Poisson2_1D.md),
[Poisson2_1D](https://inlabru-org.github.io/inlabru/reference/Poisson2_1D.md)
and
[toygroups](https://inlabru-org.github.io/inlabru/reference/toygroups.md).

## See also

Useful links:

- <http://www.inlabru.org>

- <https://inlabru-org.github.io/inlabru/>

- <https://github.com/inlabru-org/inlabru>

- Report bugs at <https://github.com/inlabru-org/inlabru/issues>

## Author

Fabian E. Bachl <bachlfab@gmail.com> and Finn Lindgren
<finn.lindgren@gmail.com>
