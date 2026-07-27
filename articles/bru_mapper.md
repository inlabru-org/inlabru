# Devel: Customised model components with the bru_mapper system

(Vignette under construction!)

## Mapper system introduction

Each inlabru latent model component generates an *effect*, given the
latent *state* vector. The purpose of the mapper system is to define the
link from state vectors to effect vectors. For an ordinary model
component, named c, this link can be represented as a matrix-vector
product \eta_c(u_c) = A_c(\text{input}) u_c, where u_c is the latent
state vector, \eta_c(u_c) is the resulting component effect vector, and
A_c(\text{input}) is a *component model matrix*. This matrix depends on
the component inputs (covariates, index information, etc) from `main`,
`group`, `replicate`, and `weights` in the component definition.

The wider scope of component definitions is discussed in the [component
vignette](https://inlabru-org.github.io/inlabru/articles/component.html).
Here, we will focus on the mapper system itself, how the basic building
blocks are used to construct the built-in mappers, and finally how to
construct new mappers.

The in-package documentation for the mapper methods is contained in four
parts:

``` r

?bru_mapper # Mapper constructors, with links to predefined mappers
?bru_mapper_generics # Generic and default methods
?bru_get_mapper # Mapper extraction methods
```

Regular users normally at most need some of the methods from
[`?bru_mapper`](https://inlabru-org.github.io/inlabru/reference/bru_mapper.md)
and links to predefined mappers, and sometimes
[`?bru_mapper_generics`](https://inlabru-org.github.io/inlabru/reference/bru_mapper_generics.md).
The
[`bru_mapper_define()`](https://inlabru-org.github.io/inlabru/reference/bru_mapper.md)
and
[`bru_get_mapper()`](https://inlabru-org.github.io/inlabru/reference/bru_get_mapper.md)
methods are needed for those implementing their own mapper class.

## Mappers

The main purpose of each mapper class is to allow evaluating a component
effect, from given `input` and latent `state`, by calling

``` r

ibm_eval(mapper, input, state)
```

for the `mapper` associated with the component definition.

### Basic mappers

Basic mappers take covariate vectors or matrices as input and numeric
vectors as state.

Constructors:

- [`bm_const()`](https://inlabru-org.github.io/inlabru/reference/bm_const.md)
- [`bm_linear()`](https://inlabru-org.github.io/inlabru/reference/bm_linear.md)
- `bm_index(n)`
- `bm_factor(values, factor_mapping)`
- [`bm_matrix()`](https://inlabru-org.github.io/inlabru/reference/bm_matrix.md)
- `bm_harmonics(order, scaling, intercept, interval)`

#### const

The `const` mapper defines a mapping from an empty state vector. This is
used to define offset components, \eta(u)\_j = z_j, where z is a fixed
input vector, indexed by j.

#### linear

The `linear` mapper defines a mapper from a length 1 state vector. This
is used to define ordinary “fixed” covariate effects, linear in both the
covariate and in the state variable; \eta(u)\_j = z_j u, where u is the
state, and z is the input covariate vector.

#### index

The `index(n)` mapper defines a direct mapping from a length `n` state
variable. This is used for structured and unstructured “random effect”
component effects; \eta(u)\_j = u\_{z_j}, where u is the state vector,
and z is the input index vector.

#### factor

The `factor(values, factor_mapping)` mapper defines a direct mapping
between a state vector of length equal to the number of factor levels of
`values` (for `factor_mapping = "full"`), or one less (for
`factor_mapping = "contrast"`). Each factor level is represented as a
dummy 0/1 variable, or equivalently, the input is used as a label index
into the state vector; \eta(u)\_j = u\_{z_j}. In the `"contrast"` case,
the first factor level is removed from the model, and the corresponding
effect defined to be zero.

#### matrix

The `matrix` mapper defines a direct matrix multiplication between a
pre-computed model matrix A\_\text{input} with the state vector; \eta(u)
= A\_\text{input} u. This is used e.g. for linear model formula input,
that is converted to a component model matrix before handing it over to
the mapper system.

#### harmonics

The `harmonics(order, scaling, intercept, interval)` mapper defines a
sum of weighted harmonics on a real interval (L, U), each frequency
additionally scaled by factors determined by `scale`; \eta(u_j) =
\sum\_{k=1}^{1+2p} A\_{j,k} u_k, where A\_{j,k}=s_k for k=1,
A\_{j,2k}=s\_{k+1}\cos\[2\pi k (z_j - L)/(U-L)\] for k=1,2,\dots,p, and
A\_{j,2k+1}=s\_{k+1}\sin\[2\pi k (z_j - L)/(U-L)\] for k=1,2,\dots,p,
where p=`order`. When `intercept` is `TRUE`, the state vector has length
1+2p. When `intercept` is `FALSE`, the first, constant, column is
removed, and the state vector has length 2p.

### Transformation mappers

Transformation mappers are mappers that would normally be combined with
other mappers, as steps in a sequence of transformations, or as
individual transformation mappers.

- [`bm_shift()`](https://inlabru-org.github.io/inlabru/reference/bm_shift.md)
- [`bm_scale()`](https://inlabru-org.github.io/inlabru/reference/bm_scale.md)
- `bm_marginal(qfun, pfun, ...)`
- `bm_aggregate(rescale)`
- `bm_logsumexp(rescale)`
- [`bm_logitaverage()`](https://inlabru-org.github.io/inlabru/reference/bm_logitaverage.md)

#### shift

The `shift` mapper adds the `input` vector to the `state` vector. This
can be used as a stage of a `pipe` mapper.

``` r

mapper <- bm_shift()
ibm_eval(mapper,
  input = ...,
  state = ...
)
```

#### scale

The `scale` mapper multiplies the `input` vector with the `state`
vector. This is used for INLA models to implement the `weights`
component scaling, as the final stage of a `pipe` mapper, see below.

``` r

mapper <- bm_scale()
ibm_eval(mapper,
  input = ...,
  state = ...
)
```

#### marginal

The `marginal` mapper transforms the `state` vector from standardised
Gaussian marginal distribution to the distribution defined by a quantile
function. This can be used to define components with latent marginal
`N(0,1)` distribution but non-Gaussian effect.

``` r

mapper <- bm_marginal(
  qfun = ..., pfun = ..., dfun = ..., ..., inverse = ...
)
ibm_eval(mapper,
  input = NULL, # If a list is given, it overrides the parameter specification
  state = ...
)

# Examples:
mapper <- bm_marginal(qfun = qexp, rate = 1 / 8)
mapper <- bm_marginal(
  qfun = qexp, pfun = pexp, dfun = dexp, rate = 1 / 8
)
```

#### aggregate

The `aggregate` mapper aggregates the `state` vector elements by
blockwise summation after scaling the elements by weights. This can be
used for summation, integration, and averaging (for `rescale=TRUE`).

``` r

mapper <- bm_aggregate(rescale = ...)
ibm_eval(mapper,
  input = list(block = ..., weights = ...),
  state = ...
)
```

#### logsumexp

The `logsumexp` mapper aggregates the `exp(state)` vector elements by
blockwise summation after scaling the elements by weights, and takes the
logarithm. The implementation takes care to avoid numerical overflow.
This can be used for summation, integration, and averaging (for
`rescale=TRUE`).

``` r

mapper <- bm_logsumexp(rescale = ...)
ibm_eval(mapper,
  input = list(block = ..., weights = ...),
  state = ...
)
```

#### logitaverage

The `logitaverage` mapper aggregates `plogis(state)` vector elements by
blockwise summation after scaling the elements by weights, and takes
[`qlogis()`](https://rdrr.io/r/stats/Logistic.html) of the results. The
implementation takes care to avoid numerical overflow. This can be used
for averaging probabilities on the logit scale. However, note that while
aggregation of Poisson variables gives a new Poisson variable,
aggregation of Binomial variables with different probability parameters
does not give a new Binomial variable, so this type of calculation
should be used with care.

``` r

mapper <- bm_logitaverage()
ibm_eval(mapper,
  input = list(block = ..., weights = ...),
  state = ...
)
```

### Compound mappers

Compound mappers define collections or chains of mappings, and can take
various forms of input. The state vector is normally a numeric vector,
but can in some cases be a list of vectors.

- `bm_collect(mappers, hidden)`
- `bm_multi(mappers)`
- `bm_pipe(mappers)`
- `bm_sum(mappers, single_input)`
- `bm_repeat(mapper, n_rep, interleaved)`

#### collect

The `collect` mapper defines a compound mapper where the Jacobian is
block-diagonal, and each block is defined by a separate sub-mapper. When
`hidden = TRUE`, it can be used for INLA models where the user-visible
part of the latent state is only a part of the internal state, such as
for `"bym"` models.

``` r

mapper <- bm_collect(list(name1 = ..., name2 = ..., ...),
  hidden = FALSE
)
ibm_eval(mapper,
  input = list(name1 = ..., name2 = ..., ...),
  state = ...
)
# If hidden = TRUE, the inla_f=TRUE argument "hides" all but the first mapper:
ibm_eval(mapper,
  input = name1_input,
  state = ...,
  inla_f = TRUE
)
```

#### multi

The `multi` mapper defines a compound mapper where the Jacobian is the
row-wise Kronecker product between the sub-mapper Jacobians. This can be
used for INLA models to construct the internal mapper for `main`,
`group`, and `replicate`, but also for user-defined models where a
single `rgeneric` or `cgeneric` model is defined on e.g. a space-time
product domain.

``` r

mapper <- bm_multi(list(name1 = ..., name2 = ..., ...))
ibm_eval(mapper,
  input = list(name1 = ..., name2 = ..., ...),
  state = ...
)
```

#### pipe

Pipe mappers chain multiple mappers into a sequence, where the
evaluation result of each mapper is given as the state vector for the
next mapper.

``` r

mapper <- bm_pipe(list(name1 = ..., name2 = ..., ...))
ibm_eval(mapper,
  input = list(name1 = ..., name2 = ..., ...),
  state = ...
)
```

#### sum (from version `2.12.0.9001`)

Sum mappers take a list of mappers and adds their effects, based on a
state vector with the individual mapper states combined in sequence. The
`ibm_` method inputs are passed on to the sub-mapper methods. When
`single_input` is `TRUE`, the same input is passed to all sub-mappers.
Otherwise, the `input` must be a `list`, `data.frame`, or `matrix`,
whose elements/columns are passed on to the respective sub-mappers.

``` r

mapper <- bm_sum(mappers = ..., single_input = ...)
ibm_eval(mapper,
  input = ...,
  state = ...
)
```

#### repeat

Repeat mappers take a single mapper and repeats it, corresponding to a
state vector of length `sum(n_rep) * ibm_n(mapper, ...)`. The `ibm_`
method inputs are passed on to the sub-mapper methods.

If `interleaved` is `TRUE` (from version `2.12.0.9001`), the repeated
mappers are interleaved, e.g. if the state vectors for two individual
mappers were `x1` and `x2`, interleaved mappers have state vector
`c(x1[1], x2[1], ..., x1[2], x2[2], ...)`, whereas ordinary
non-interleaved mappers have state vector `c(x1, x2, ...)`.

For version `2.12.0.9001`, `n_rep` and `interleaved` may be vectors, for
grouped repetition and interleaving. Non-interleaved repeats are
combined, so that the resulting internal mapper structure is a compact
as possible.

``` r

mapper <- bm_repeat(mapper = ..., n_rep = ..., interleaved = ...)
ibm_eval(mapper,
  input = ...,
  state = ...
)
```

#### The core model component mapper

All model component mappers are currently defined as a `pipe` mapper
containing a `multi` mapper followed by a `scale` mapper:

``` r

mapper <-
  bm_pipe(
    list(
      mapper = bm_multi(list(main = ..., group = ..., replicate = ...)),
      scale = bm_scale()
    )
  )
ibm_eval(mapper,
  input = list(
    mapper = list(main = ..., group = ..., replicate = ...),
    scale = ...
  ),
  state = ...
)
```

### Object mappers

For object with predefined mapper conversion methods, use
`bru_mapper(object)` to obtain a suitable mapper object. Current
predefined object mappers are defined for objects of these classes:

- `fm_mesh_2d`; Use `bru_mapper(object)` or `bm_fmesher(object)`
- `fm_mesh_1d`; Use `bru_mapper(object, indexed)` or, if `indexed=TRUE`,
  `bm_fmesher(object)`.
- Any object supporting the
  [`fmesher::fm_dof()`](https://inlabru-org.github.io/fmesher/reference/fm_dof.html)
  and
  [`fmesher::fm_basis()`](https://inlabru-org.github.io/fmesher/reference/fm_basis.html)
  methods; Use `bm_fmesher(object)`. If `fmesher` in the future gains a
  common base class for such objects, a corresponding
  [`bru_mapper()`](https://inlabru-org.github.io/inlabru/reference/bru_mapper.md)
  conversion method will be added.

More details are given below.

### Model object mappers

For inla model objects with predefined mappers, use
`bru_get_mapper(object)` Current predefined model object mappers are
defined for models of these classes:

- `inla.spde`; this includes
  [`inla.spde2.matern()`](https://rdrr.io/pkg/INLA/man/inla.spde2.matern.html)
  and
  [`inla.spde2.pcmatern()`](https://rdrr.io/pkg/INLA/man/inla.spde2.pcmatern.html)
  models. The conversion method calls
  [`bru_mapper()`](https://inlabru-org.github.io/inlabru/reference/bru_mapper.md)
  in the appropriate way for the mesh type used for the model
- `inla.rgeneric`; This relies on the rgeneric definition including a
  `"mapper"` callback argument. A less invasive method, that works for
  both rgeneric and cgeneric models, is to define a
  [`bru_get_mapper()`](https://inlabru-org.github.io/inlabru/reference/bru_get_mapper.md)
  method for the class, which needs to include a unique class
  identifier, e.g. add
  `class(object) <- c("my_unique_model_class", class(object))` and
  define

``` r

bru_get_mapper.my_unique_model_class <- function(model, ...) {
  ...
}
```

and register it in the namespace.

### Special mappers

- [`bm_taylor()`](https://inlabru-org.github.io/inlabru/reference/bm_taylor.md)

## Mapper methods

``` r

ibm_n(mapper, inla_f, ...)
ibm_n_output(mapper, input, ...)
ibm_values(mapper, inla_f, ...)
ibm_jacobian(mapper, input, state, ...)
ibm_eval(mapper, input, state, ...)
ibm_names(mapper, ...)
ibm_inla_subset(mapper, ...)
```

An example where the `inla_f` argument matters is the `bm_collect`
class, when the `hidden=TRUE` argument is used to indicate that only the
first mapper should be used for the
[`INLA::f()`](https://rdrr.io/pkg/INLA/man/f.html) inputs, e.g. for
`"bym2"` models. For `inla_f=FALSE` (the default), the `ibm_n` and
`ibm_values` methods return the total number of latent variables, but
with `inla_f=TRUE` we get the values needed by
[`INLA::f()`](https://rdrr.io/pkg/INLA/man/f.html) instead:

``` r

mapper <- bm_collect(
  list(
    a = bm_index(3),
    b = bm_index(2)
  ),
  hidden = TRUE
)
ibm_n(mapper)
#> [1] 5
ibm_values(mapper)
#> [1] 1 2 3 4 5
ibm_n(mapper, inla_f = TRUE)
#> [1] 3
ibm_values(mapper, inla_f = TRUE)
#> [1] 1 2 3
```

``` r

ibm_n(mapper, multi = TRUE)
#> $a
#> [1] 3
#> 
#> $b
#> [1] 2
ibm_values(mapper, multi = TRUE)
#> $a
#> [1] 1 2 3
#> 
#> $b
#> [1] 1 2
ibm_n(mapper, inla_f = TRUE, multi = TRUE)
#> $a
#> [1] 3
#> 
#> $b
#> [1] 2
ibm_values(mapper, inla_f = TRUE, multi = TRUE)
#> $a
#> [1] 1 2 3
#> 
#> $b
#> [1] 1 2
```

For the `bm_multi` class, the `multi` argument provides access to the
inner layers of the multi-mapper:

``` r

mapper <- bm_multi(list(
  a = bm_index(3),
  b = bm_index(2)
))
ibm_n(mapper)
#> [1] 6
ibm_n(mapper, multi = TRUE)
#> $a
#> [1] 3
#> 
#> $b
#> [1] 2
ibm_values(mapper)
#> [1] 1 2 3 4 5 6
ibm_values(mapper, multi = TRUE)
#>   a b
#> 1 1 1
#> 2 2 1
#> 3 3 1
#> 4 1 2
#> 5 2 2
#> 6 3 2
```

``` r

ibm_n(mapper, inla_f = TRUE)
#> [1] 6
ibm_n(mapper, multi = TRUE, inla_f = TRUE)
#> $a
#> [1] 3
#> 
#> $b
#> [1] 2
ibm_values(mapper, inla_f = TRUE)
#> [1] 1 2 3 4 5 6
ibm_values(mapper, multi = TRUE, inla_f = TRUE)
#>   a b
#> 1 1 1
#> 2 2 1
#> 3 3 1
#> 4 1 2
#> 5 2 2
#> 6 3 2
```

The default `ibm_inla_subset` method determines the inla subset by
comparing the full values information from
`ibm_values(mapper, inla_f = FALSE)` to the inla specific values
information from `ibm_values(mapper, inla_f = TRUE)`, to determine the
logical vector identifying the inla values subset. Custom mappers should
normally not need to specialise this method.

## Mappers for `fmesher` objects

For component models referenced by a character label (e.g. `"iid"`,
`"bym2"`, `"rw2"` etc), inlabru will construct default mappers, that in
most cases replicate the default
[`INLA::f()`](https://rdrr.io/pkg/INLA/man/f.html) behaviour.

All `fmesher` objects with defined
[`fmesher::fm_dof()`](https://inlabru-org.github.io/fmesher/reference/fm_dof.html)
and
[`fmesher::fm_basis()`](https://inlabru-org.github.io/fmesher/reference/fm_basis.html)
methods are supported by the
[`bm_fmesher()`](https://inlabru-org.github.io/inlabru/reference/bm_fmesher.md)
constructor:

``` r

bm_fmesher(mesh)
```

This defines a mapper that gives `fm_dof(mesh)` as output from
[`ibm_n()`](https://inlabru-org.github.io/inlabru/reference/ibm_n.md),
and `seq_len(fm_dof(mesh))` as output from
[`ibm_values()`](https://inlabru-org.github.io/inlabru/reference/ibm_values.md)
(a so-called “indexed mapper”).

For the `fm_mesh_1d` and `fm_mesh_2d` classes, the default mappers can
also be constructed by `bru_mapper` method:

``` r

bru_mapper(mesh)
```

For 1d meshes (`fm_mesh_1d`), one can specify a non-indexed mapper with
the `indexed = FALSE` argument, that gives `mesh$mid` or `mesh$loc` as
output

``` r

# If ibm_values() should return
#   mesh$mid (if available, e.g. for "rw2" models with degree=2 meshes)
# or
#   mesh$loc
bru_mapper(mesh, indexed = FALSE)
# If ibm_values() should return seq_len(fm_dof(mesh))
#   (e.g. for inla.spde2.pcmatern() models)
bru_mapper(mesh, indexed = TRUE)
```

## Customised mappers

A mapper object should store enough information in order for the `ibm_*`
methods to work. The simplest case of a customised mapper is to just
attached a new class label to the front of the S3
[`class()`](https://rdrr.io/r/base/class.html) information of an
existing mapper, to obtain a class to override some of the standard
`ibm_*` method implementations. More commonly, one would instead start
with a basic [`list()`](https://rdrr.io/r/base/list.html), that might
*contain* an existing mapper object, and then add methods that know how
to use that information. In all these cases, the
[`bru_mapper_define()`](https://inlabru-org.github.io/inlabru/reference/bru_mapper.md)
method should be used to properly set the class information:

``` r

bru_mapper_define(mapper, new_class)
```

Implementations can avoid having to define `ibm_n` and `ibm_values`
methods, by instead computing and storing `n`, `n_inla`, `values`, and
`values_inla` in the mapper object during construction:

``` r

bru_mapper_define(
  mapper = list(n = 10, values = 1:10),
  new_class = "my_bm_class_name"
)
```

The default `ibm_n` and `ibm_values` methods checks if these values are
available, and return the appropriate values depending on the `inla_f`
argument. When `*_inla` values are requested but not available, the
methods fall back to the non-inla versions. If the needed information is
not found, the default methods give an error message.

Note: Before version 2.6.0, `ibm_amatrix` was used instead of
`ibm_jacobian`. Version 2.6.0 still supports this, by having the default
`ibm_jacobian` method call `ibm_amatrix`, but will give a deprecation
message later, and it may be removed in a future version.

## Example

Let’s build a mapper class `bm_p_quadratic` for a model component that
takes input covariates with values a\_{ij}, i=1,\dots,m and j=1,\dots,p,
in a `matrix` or `data.frame` and evaluates a full quadratic expression
\eta_i = x_0 + \sum\_{j=1}^p a\_{ij} x_j +
\frac{1}{2}\sum\_{j=1}^p\sum\_{k=1}^j \gamma\_{j,k} a\_{ij}a\_{ik}
x\_{j,k}, where the latent component vector is
\boldsymbol{x}=\[x_0,x_1,\dots,x_p,x\_{1,1},\dots,x\_{p,1},x\_{2,2}\dots,x\_{p,p}\],
and \gamma\_{j,k} = \begin{cases} 1, & j = k,\\ 2, & j \> k. \end{cases}

We start with the constructor method. Like the `bm_matrix`, we require
the user to supply a vector of covariate labels, and store them as a
character vector. We also include a `min_degree` parameter to control if
an intercept (`min_degree <= 0`) and linear terms (`min_degree <= 1`)
should be included in the model.

``` r

bm_p_quadratic <- function(labels, min_degree = 0, ...) {
  if (is.factor(labels)) {
    mapper <- list(
      labels = levels(labels),
      min_degree = min_degree
    )
  } else {
    mapper <- list(
      labels = as.character(labels),
      min_degree = min_degree
    )
  }
  bru_mapper_define(mapper, new_class = "bm_p_quadratic")
}
```

The `ibm_n` method can compute the value of n from
n=1+p+\frac{p(p+1)}{2}:

``` r

ibm_n.bm_p_quadratic <- function(mapper, ...) {
  p <- length(mapper$labels)
  (mapper$min_degree <= 0) + (mapper$min_degree <= 1) * p + p * (p + 1) / 2
}
```

For `ibm_values`, the default method is sufficient, returning the vector
\[1,\dots,n\] with n obtained by `ibm_n(mapper)`. However, for clearer
result naming, we can use a character vector, at least for INLA
[`f()`](https://rdrr.io/pkg/INLA/man/f.html) models that allow it (we
could have an option argument to the `bm_p_quadratic` constructor to
control this):

``` r

ibm_values.bm_p_quadratic <- function(mapper, ...) {
  p <- length(mapper$labels)
  n <- ibm_n(mapper)
  jk <- expand.grid(seq_len(p), seq_len(p))
  jk <- jk[jk[, 2] <= jk[, 1], , drop = FALSE]
  c(
    if (mapper$min_degree <= 0) "Intercept" else NULL,
    if (mapper$min_degree <= 1) mapper$labels else NULL,
    paste0(mapper$labels[jk[, 1]], ":", mapper$labels[jk[, 2]])
  )
}
```

While the approach to `ibm_n` and `ibm_values` above works, it is
wasteful to recompute the `n` and `values` information for each call.
Instead we can use that the default method will check if `n` or `values`
are stored in the mapper object, and then return them. If we change the
`.` in the method definitions to `_`, we can end the constructor method
like this, making subsequent method calls faster:

``` r

bm_p_quadratic <- function(labels, min_degree = 0, ...) {
  ...
  mapper <- bru_mapper_define(mapper, new_class = "bm_p_quadratic")
  mapper$n <- ibm_n_bm_p_quadratic(mapper)
  mapper$values <- ibm_values_bm_p_quadratic(mapper)
  mapper
}
```

Note that the order matters, since our `ibm_values_bm_p_quadratic`
function calls `ibm_n(mapper)`.

We can now define the main part of the mapper interface that computes
the model matrix linking the latent variables to the component effect.
It’s required that `NULL` input should return a 0-by-n matrix.

``` r

ibm_jacobian.bm_p_quadratic <- function(mapper, input, ...) {
  if (is.null(input)) {
    return(Matrix::Matrix(0, 0, ibm_n(mapper)))
  }
  p <- length(mapper$labels)
  n <- ibm_n(mapper)
  N <- NROW(input)
  A <- list()
  in_ <- as(input, "Matrix")
  idx <- 0
  if (mapper$min_degree <= 0) {
    idx <- idx + 1
    A[[idx]] <- Matrix::Matrix(1, N)
  }
  if (mapper$min_degree <= 1) {
    idx <- idx + 1
    A[[idx]] <- in_
  }
  for (k in seq_len(p)) {
    idx <- idx + 1
    A[[idx]] <- in_[, seq(k, p, by = 1), drop = FALSE] * in_[, k]
    A[[idx]][, k] <- A[[idx]][, k] / 2
  }
  A <- do.call(cbind, A)
  colnames(A) <- as.character(ibm_values(mapper))
  A
}
```

If the output of a mapper is of a different size than `NROW(input)`, the
mapper should define a `ibm_n_output(mapper, input, ...)` method to
return the output size for a given input.
