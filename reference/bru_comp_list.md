# Methods for inlabru component lists

Constructor methods for inlabru component lists. Syntax details are
given in
[`bru_comp()`](https://inlabru-org.github.io/inlabru/reference/bru_comp.md).

## Usage

``` r
bru_comp_list(object, ..., .envir = parent.frame())

# S3 method for class 'formula'
bru_comp_list(object, ..., lhoods = NULL, .envir = parent.frame())

# S3 method for class 'list'
bru_comp_list(
  object,
  ...,
  lhoods = NULL,
  .envir = parent.frame(),
  inputs = NULL
)

# S3 method for class 'bru_comp'
bru_comp_list(object, ..., .envir = parent.frame())

# S3 method for class 'bru_comp_list'
bru_comp_list(object, ..., .envir = parent.frame())

# S3 method for class 'bru_comp_list'
c(...)

# S3 method for class 'bru_comp'
c(...)

# S3 method for class 'bru_comp_list'
x[i]
```

## Arguments

- object:

  The object to operate on

- ...:

  Parameters passed on to other methods. Also see Details.

- .envir:

  An evaluation environment for non-formula input

- lhoods:

  A
  [bru_obs_list](https://inlabru-org.github.io/inlabru/reference/bru_obs.md)
  object

- inputs:

  A tree-like list of component input evaluations, from
  [`bru_input.bru_obs_list()`](https://inlabru-org.github.io/inlabru/reference/bru_input.md).

- x:

  `bru_comp_list` object from which to extract a sub-list

- i:

  indices specifying elements to extract

## Value

A `bru_comp_list` object, which is a named list of `bru_comp` objects.

## Methods (by class)

- `bru_comp_list(formula)`: Convert a component formula into a
  `bru_comp_list` object

- `bru_comp_list(list)`: Combine a list of components, component lists,
  and/or component formulas into a single `bru_comp_list` object

- `bru_comp_list(bru_comp)`: Place a single `bru_comp` object into a
  `bru_comp_list` object.

- `bru_comp_list(bru_comp_list)`: Make sure a `bru_comp_list` object is
  fully configured.

## Methods (by generic)

- `c(bru_comp_list)`: The `...` arguments should be `bru_comp_list`
  objects. The environment from the first argument will be applied to
  the resulting `bru_comp_list`.

## Functions

- `c(bru_comp)`: The `...` arguments should be `component` objects from
  [`bru_comp()`](https://inlabru-org.github.io/inlabru/reference/bru_comp.md).
  The environment from the first argument will be applied to the
  resulting `bru_comp_list`.

## See also

Other component constructors:
[`bru_comp()`](https://inlabru-org.github.io/inlabru/reference/bru_comp.md)

## Author

Fabian E. Bachl <bachlfab@gmail.com> and Finn Lindgren
<finn.lindgren@gmail.com>

## Examples

``` r
# As an example, let us create a linear component. Here, the component is
# called "myLinearEffectOfX" while the covariate the component acts on is
# called "x". Note that a list of components is returned because the
# formula may define multiple components

eff <- bru_comp_list(~ myLinearEffectOfX(main = x, model = "linear"))
summary(eff[[1]])
#> Label:   myLinearEffectOfX 
#>   Type:  main = linear 
#>   Map:   pipe = multi(main = {autodetect(x)}) 
#>   INLA formula:  
#>     ~ . + f(myLinearEffectOfX, model =
#>       BRU_myLinearEffectOfX_main_model) 
# Equivalent shortcuts:
eff <- bru_comp_list(~ myLinearEffectOfX(x, model = "linear"))
eff <- bru_comp_list(~ myLinearEffectOfX(x))
# Individual component
eff <- bru_comp("myLinearEffectOfX", main = x, model = "linear")
```
