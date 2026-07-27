# Constructors for `bru_mapper` objects

Constructors for `bru_mapper` objects

## Usage

``` r
bru_mapper(...)

bru_mapper_define(mapper, new_class = NULL, remove_class = "list", ...)
```

## Arguments

- ...:

  Arguments passed on to sub-methods, or used for special purposes, see
  details for each function below.

- mapper:

  For `bru_mapper_define`, a prototype mapper object, see Details.

- new_class:

  If non-`NULL`, this is added at the front of the class definition

- remove_class:

  If non-`NULL`, this class or classes is removed from the class
  definition before adding the `new_class` names. Default is `"list"`.

## Value

A `bru_mapper` object.

## Functions

- `bru_mapper()`: Generic mapper S3 constructor, used for constructing
  mappers for special objects. See below for details of the default
  constructor `bru_mapper_define()` that can be used to define new
  mapper clesses in user code. To extracting mappers for latent
  component models, see
  [`bru_get_mapper()`](https://inlabru-org.github.io/inlabru/reference/bru_get_mapper.md).

- `bru_mapper_define()`: Adds the `new_class` and `"bru_mapper"` class
  names to the inheritance list for the input `mapper` object, unless
  the object already inherits from these.

  To register mapper classes and methods in scripts, use
  [`.S3method()`](https://rdrr.io/r/base/S3method.html) to register the
  methods, e.g.
  `.S3method("ibm_jacobian", "my_mapper_class", ibm_jacobian.my_mapper_class)`.

  In packages with `Suggests: inlabru`, add method information for
  delayed registration, e.g.:

      #' @rawNamespace S3method(inlabru::bru_get_mapper, inla_rspde)
      #' @rawNamespace S3method(inlabru::ibm_n, bru_mapper_inla_rspde)
      #' @rawNamespace S3method(inlabru::ibm_values, bru_mapper_inla_rspde)
      #' @rawNamespace S3method(inlabru::ibm_jacobian, bru_mapper_inla_rspde)

  or before each method, use `@exportS3Method`:

      #' @exportS3Method inlabru::bru_get_mapper

  etc., which semi-automates it.

## See also

[bru_mapper_generics](https://inlabru-org.github.io/inlabru/reference/bru_mapper_generics.md)
for generic methods, the individual mapper pages for special method
implementations, and
[bru_get_mapper](https://inlabru-org.github.io/inlabru/reference/bru_get_mapper.md)
for hooks to extract mappers from latent model object class objects.

Other mappers:
[`bm_aggregate()`](https://inlabru-org.github.io/inlabru/reference/bm_aggregate.md),
[`bm_collect()`](https://inlabru-org.github.io/inlabru/reference/bm_collect.md),
[`bm_const()`](https://inlabru-org.github.io/inlabru/reference/bm_const.md),
[`bm_expr()`](https://inlabru-org.github.io/inlabru/reference/bm_expr.md),
[`bm_factor()`](https://inlabru-org.github.io/inlabru/reference/bm_factor.md),
[`bm_fm_mesh_1d`](https://inlabru-org.github.io/inlabru/reference/bm_fm_mesh_1d.md),
[`bm_fmesher()`](https://inlabru-org.github.io/inlabru/reference/bm_fmesher.md),
[`bm_harmonics()`](https://inlabru-org.github.io/inlabru/reference/bm_harmonics.md),
[`bm_index()`](https://inlabru-org.github.io/inlabru/reference/bm_index.md),
[`bm_linear()`](https://inlabru-org.github.io/inlabru/reference/bm_linear.md),
[`bm_logitaverage()`](https://inlabru-org.github.io/inlabru/reference/bm_logitaverage.md),
[`bm_logsumexp()`](https://inlabru-org.github.io/inlabru/reference/bm_logsumexp.md),
[`bm_marginal()`](https://inlabru-org.github.io/inlabru/reference/bm_marginal.md),
[`bm_matrix()`](https://inlabru-org.github.io/inlabru/reference/bm_matrix.md),
[`bm_multi()`](https://inlabru-org.github.io/inlabru/reference/bm_multi.md),
[`bm_pipe()`](https://inlabru-org.github.io/inlabru/reference/bm_pipe.md),
[`bm_reparam()`](https://inlabru-org.github.io/inlabru/reference/bm_reparam.md),
[`bm_repeat()`](https://inlabru-org.github.io/inlabru/reference/bm_repeat.md),
[`bm_scale()`](https://inlabru-org.github.io/inlabru/reference/bm_scale.md),
[`bm_shift()`](https://inlabru-org.github.io/inlabru/reference/bm_shift.md),
[`bm_sum()`](https://inlabru-org.github.io/inlabru/reference/bm_sum.md),
[`bm_taylor()`](https://inlabru-org.github.io/inlabru/reference/bm_taylor.md),
[`bru_get_mapper()`](https://inlabru-org.github.io/inlabru/reference/bru_get_mapper.md)

## Examples

``` r
mapper <- bm_index(5)
ibm_jacobian(mapper, input = c(1, 3, 4, 5, 2))
#> 5 x 5 sparse Matrix of class "dgCMatrix"
#>               
#> [1,] 1 . . . .
#> [2,] . . 1 . .
#> [3,] . . . 1 .
#> [4,] . . . . 1
#> [5,] . 1 . . .
```
