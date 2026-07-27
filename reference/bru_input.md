# Store and evaluate component inputs

Store and evaluate component inputs, including support for
non-standard-evaluation with `rlang`, automatic transfer to function
calls,
[`eval_spatial()`](https://inlabru-org.github.io/inlabru/reference/eval_spatial.md),
and
[`MatrixModels::model.Matrix()`](https://rdrr.io/pkg/MatrixModels/man/model.Matrix.html).

These functions are normally only called internally, but users may find
it useful to call them directly to experiment, and to make advanced
model definitions using
[ibm_input_set](https://inlabru-org.github.io/inlabru/reference/ibm_input.md)/[`ibm_input_new()`](https://inlabru-org.github.io/inlabru/reference/ibm_input.md).

## Usage

``` r
new_bru_input(input, label = NULL, layer = NULL, selector = NULL, ...)

bru_input(...)

# S3 method for class 'bru_input'
bru_input(
  x,
  data = NULL,
  mask = NULL,
  null.on.fail = FALSE,
  .envir = parent.frame(),
  ...
)

# S3 method for class 'bru_comp'
bru_input(x, ..., label = x$label)

# S3 method for class 'bru_comp_list'
bru_input(x, ...)

# S3 method for class 'bru_model'
bru_input(x, lhoods = deprecated(), ...)

# S3 method for class 'bru_obs'
bru_input(x, components, ...)

# S3 method for class 'bru_obs_list'
bru_input(x, components, ...)

# S3 method for class 'bru_mapper'
bru_input(x, ..., label = "<unknown>")

# S3 method for class 'bm_pipe'
bru_input(x, ..., label = "<unknown>")

# S3 method for class 'bm_multi'
bru_input(x, ..., label = "<unknown>")

# S3 method for class 'bm_collect'
bru_input(x, ..., label = "<unknown>")

# S3 method for class 'bm_repeat'
bru_input(x, ..., label = "<unknown>")

# S3 method for class 'bm_sum'
bru_input(x, ..., label = "<unknown>")
```

## Arguments

- input:

  An expression to be evaluated.

- label:

  character; optional label used to identify the object in informational
  messages.

- layer:

  Optional expression that evaluates layer names or indices for use with
  [`eval_spatial()`](https://inlabru-org.github.io/inlabru/reference/eval_spatial.md)

- selector:

  character or integer; optional selector for use with use with
  [`eval_spatial()`](https://inlabru-org.github.io/inlabru/reference/eval_spatial.md)

- ...:

  Passed on to sub-methods.

- x:

  A `bru_input` object, or other object for recursive evaluation.

- data:

  A `data.frame`, `tibble`, `sf`, `list`, or `Spatial*` object of
  covariates and/or point locations.

- mask:

  A `bru_data_mask` object, constructed internally.

- null.on.fail:

  logical; if `TRUE`, return `NULL` if the input evaluation fails. If
  `FALSE` (default), return a vector of 1s and issue a warning (only for
  deprecated use of `coordinates` without having loaded the `sp`
  package), or stop with an error.

- .envir:

  environment in which to evaluate the input expression. Default is
  [`parent.frame()`](https://rdrr.io/r/base/sys.parent.html)

- lhoods:

  **\[deprecated\]** from `2.14.1.9011`, as it's now part of the
  [bru_model](https://inlabru-org.github.io/inlabru/reference/bru_model.md)
  object. A `bru_obs_list` object containing all observations models.

- components:

  A `bru_comp_list` object containing all components defined in the
  model.

## Value

- `bru_input(bru_comp)`: A list of mapper input values, formatted for
  the full component mapper (of type
  [bm_pipe](https://inlabru-org.github.io/inlabru/reference/bm_pipe.md))

&nbsp;

- `bru_input(bru_comp_list)`: A list of mapper input values, with one
  entry for each component.

&nbsp;

- `bru_input(bru_model)`: A list of mapper input values, with one entry
  for each observation model, each containing a list of inputs for the
  components used by the corresponding observation model.

&nbsp;

- `bru_input(bru_obs)`: A list of mapper input values, for each of the
  components used by the corresponding observation model.

&nbsp;

- `bru_input(bru_obs_list)`: A list of mapper input values, with one
  entry for each observation model, each containing a list of inputs for
  the components used by the corresponding observation model.

&nbsp;

- `bru_input(bru_mapper)`: The evaluated input associated with a mapper;
  see
  [ibm_input](https://inlabru-org.github.io/inlabru/reference/ibm_input.md)
  for details on how to associate a `bru_input` with a `bru_mapper`.

## Methods (by class)

- `bru_input(bru_input)`: Attempts to evaluate a component input (e.g.
  `main`, `group`, `replicate`, or `weight`), and process the results:

  1.  If
      [`rlang::eval_tidy()`](https://rlang.r-lib.org/reference/eval_tidy.html)
      failed, return NULL or map everything to 1 (see the `null.on.fail`
      argument). This should normally not happen, unless the component
      use logic is incorrect, (via `used`) leading to missing columns
      for a certain likelihood in a
      multi-[`bru_obs()`](https://inlabru-org.github.io/inlabru/reference/bru_obs.md)
      model.

  2.  If we obtain a function, apply the function to the data object

  3.  If we obtain an object supported by
      [`eval_spatial()`](https://inlabru-org.github.io/inlabru/reference/eval_spatial.md),
      extract the values of that data frame at the point locations

  4.  If we obtain a formula, call `ModelMatrix::model.Matrix()`

  5.  Else we obtain a vector and return as-is. This happens when input
      references a column of the data points, or some other complete
      expression

- `bru_input(bru_model)`: Computes the component inputs for included
  components for each observation model.

- `bru_input(bru_mapper)`: Evaluate the input associated with a
  `bru_mapper`.

- `bru_input(bm_pipe)`: Evaluate the inputs for each sub-mapper in a
  `bm_pipe` object.

- `bru_input(bm_multi)`: Evaluate the inputs for each sub-mapper in a
  `bm_multi` object.

- `bru_input(bm_collect)`: Evaluate the inputs for each sub-mapper in a
  `bm_collect` object.

- `bru_input(bm_repeat)`: Evaluate the inputs for the sub-mapper in a
  `bm_repeat` object.

- `bru_input(bm_sum)`: Evaluate the inputs for each sub-mapper in a
  `bm_sum` object.

## Functions

- `new_bru_input()`: Create a `bru_input` object.

- `bru_input()`: Evaluate inputs

## Simple covariates and the input parameters

It is not unusual for a random effect act on a transformation of a
covariate. In other frameworks this would mean that the transformed
covariate would have to be calculated in advance and added to the data
frame that is usually provided via the `data` parameter. inlabru
provides the option to do this transformation automatically. For
instance, one might be interested in the effect of a covariate \\x^2\\.
In inla and other frameworks this would require to add a column
`xsquared` to the input data frame and use the formula

- `formula = y ~ f(xsquared, model = "linear")`,

In inlabru this can be achieved in several ways of using the `main`
parameter (`map` in version 2.1.13 and earlier), which does not need to
be named.

- `components = y ~ psi(x^2, model = "linear")`

- `components = y ~ psi(main = x^2, model = "linear")`

- `components = y ~ psi(mySquareFun(x), model = "linear")`,

- `components = y ~ psi(myOtherSquareFun, model = "linear")`,

In the first example inlabru will interpret the map parameter as an
expression to be evaluated within the data provided. Since \\x\\ is a
known covariate it will know how to calculate it. The second example is
an expression as well but it uses a function called `mySquareFun`. This
function is defined by user but has to be accessible within the work
space when setting up the components. The third example provides the
function `myOtherSquareFun`. In this case, inlabru will call the
function as `myOtherSquareFun(.data.)`, where `.data.` is the data
provided via the
[`bru_obs()`](https://inlabru-org.github.io/inlabru/reference/bru_obs.md)
`data` parameter. The function needs to know what parts of the data to
use to construct the needed output. For example,

    myOtherSquareFun <- function(data) {
      data[ ,"x"]^2
    }

Interactions can be handled by a formula input and `model = "fixed"`:
`components = y ~ 0 + name(~ 1 + x:z, model = "fixed")`

## Spatial Covariates

When fitting spatial models it is common to work with covariates that
depend on space, e.g. sea surface temperature or elevation. Although it
is straightforward to add this data to the input data frame or write a
covariate function like in the previous section there is an even more
convenient way in inlabru. Spatial covariates are often stored as
`SpatialPixelsDataFrame`, `SpatialPixelsDataFrame` or `RasterLayer`
objects. These can be provided directly via the input expressions if
they are supported by
[`eval_spatial()`](https://inlabru-org.github.io/inlabru/reference/eval_spatial.md),
and the
[`bru_obs()`](https://inlabru-org.github.io/inlabru/reference/bru_obs.md)
data is an `sf` or `SpatialPointsDataFrame` object. `inlabru` will then
automatically evaluate and/or interpolate the covariate at your data
locations when using code like

    components = y ~ psi(mySpatialPixels, model = "linear")

For more precise control, use the the `layer` and `selector` arguments
(see
[`bru_comp()`](https://inlabru-org.github.io/inlabru/reference/bru_comp.md)),
or call
[`eval_spatial()`](https://inlabru-org.github.io/inlabru/reference/eval_spatial.md)
directly, e.g.:

    components = y ~ psi(eval_spatial(mySpatialPixels, where = .data.),
                         model = "linear")

## Coordinates

A common spatial modelling component when using inla are SPDE models. An
important feature of inlabru is that it will automatically calculate the
so called A-matrix (a component model matrix) which maps SPDE values at
the mesh vertices to values at the data locations. For this purpose, the
input can be set to `coordinates`, which is the `sp` package function
that extracts point coordinates from the `SpatialPointsDataFrame` that
was provided as input to
[`bru_obs()`](https://inlabru-org.github.io/inlabru/reference/bru_obs.md).
The code for this would look as follows:

    components = y ~ field(coordinates, model = inla.spde2.matern(...))

Since `coordinates` is a function from the `sp` package, this results in
evaluation of `sp::coordinates(.data.)`, which loses any CRS information
from the data object.

For `sf` data with a geometry column (by default named `geometry`), use

    components = y ~ field(geometry, model = inla.spde2.matern(...))

Since the CRS information is part of the geometry column of the `sf`
object, this retains CRS information, so this is more robust, and allows
the model to be built on a different CRS than the observation data.

## See also

[`bru_input_text()`](https://inlabru-org.github.io/inlabru/reference/bru_input_text.md)

[`summary.bru_input()`](https://inlabru-org.github.io/inlabru/reference/summary.bru_input.md),
[`bru_comp()`](https://inlabru-org.github.io/inlabru/reference/bru_comp.md)

[ibm_input](https://inlabru-org.github.io/inlabru/reference/ibm_input.md)

## Author

Fabian E. Bachl <bachlfab@gmail.com>, Finn Lindgren
<finn.lindgren@gmail.com>

## Examples

``` r
(inp <- new_bru_input(x, "LABEL"))
#> LABEL = x
bru_input(inp, data.frame(x = 1:3))
#> [1] 1 2 3
```
