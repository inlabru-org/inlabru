# Create an inlabru model object from model components

The
[inlabru](https://inlabru-org.github.io/inlabru/reference/inlabru-package.md)
syntax for model formulae is different from what
[`INLA::inla`](https://rdrr.io/pkg/INLA/man/inla.html) considers a
valid. In inla most of the effects are defined by adding an `f(...)`
expression to the formula. In
[inlabru](https://inlabru-org.github.io/inlabru/reference/inlabru-package.md)
the `f` is replaced by an arbitrary (exceptions: `const` and `offset`)
string that will determine the label of the effect. See Details for
further information.

## Usage

``` r
bru_model(components, lhoods, options = list(), .envir = parent.frame())

# S3 method for class 'bru_model'
summary(object, ...)

# S3 method for class 'summary_bru_model'
print(x, ...)

# S3 method for class 'bru_model'
print(x, ...)
```

## Arguments

- components:

  A
  [bru_comp_list](https://inlabru-org.github.io/inlabru/reference/bru_comp_list.md)
  object

- lhoods:

  Either a
  [`bru_obs()`](https://inlabru-org.github.io/inlabru/reference/bru_obs.md)
  object, a
  [`bru_obs_list()`](https://inlabru-org.github.io/inlabru/reference/bru_obs.md)
  object, or a list of one or more
  [`bru_obs()`](https://inlabru-org.github.io/inlabru/reference/bru_obs.md)
  or
  [`bru_obs_list()`](https://inlabru-org.github.io/inlabru/reference/bru_obs.md)
  objects, or a list of arguments for a call to
  [`bru_obs()`](https://inlabru-org.github.io/inlabru/reference/bru_obs.md).

- options:

  A
  [bru_options](https://inlabru-org.github.io/inlabru/reference/bru_options.md)
  options object or a list of options passed on to
  [`bru_options()`](https://inlabru-org.github.io/inlabru/reference/bru_options.md)

- .envir:

  The environment in which the components are evaluated.

- object:

  Object to operate on

- ...:

  Arguments passed on to other methods

- x:

  An object to be printed

## Value

A bru_model object

## Details

For instance

`y ~ f(myspde, ...)`

in INLA is equivalent to

`y ~ myspde(...)`

in inlabru.

A disadvantage of the inla way is that there is no clear separation
between the name of the covariate and the label of the effect.
Furthermore, for some models like SPDE it is much more natural to use
spatial coordinates as covariates rather than an index into the SPDE
vertices. For this purpose
[inlabru](https://inlabru-org.github.io/inlabru/reference/inlabru-package.md)
provides the `main` argument. For convenience, the `main` argument can
be used like the first argument of the f function, e.g., and is the
first argument of the component definition. The `INLA` model formula

`y ~ f(temperature, model = 'linear')`

is equivalent to the `inlabru` component and formula definition

`y ~ temperature(temperature, model = 'linear')` and
`y ~ temperature(main = temperature, model = 'linear')` as well as
`y ~ temperature(model = 'linear')` which sets `main = temperature`.

On the other hand, `main` can also be a function mapping, e.g the
[`sp::coordinates()`](https://edzer.github.io/sp/reference/coordinates.html)
function:

`y ~ mySPDE(coordinates, ...)`

This extracts the coordinates from the data object, and maps it to the
latent field via the information given in the `mapper`, which by default
is extracted from the `model` object, in the case of `spde` model
objects.

Morevover, `main` can be any expression that evaluates within your data
as an environment. For instance, if your data has columns 'a' and 'b',
you can create a fixed effect of 'sin(a+b)' by setting `map` in the
following way:

`y ~ myEffect(sin(a+b))`
