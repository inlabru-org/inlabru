# Create inlabru data mask

Create a data mask for inlabru, which allows for evaluating expressions
in the context of the data objects, with support for tidy evaluation
pronouns and direct access to the data objects as `.name.`. This is an
internal inlabru method, not intended for general use.

## Usage

``` r
bru_data_mask(data, pronouns = NULL, objects = NULL)
```

## Arguments

- data:

  list of data objects in priority order. Named elements will be
  available as whole objects `.name.` as well as pronouns `.name` in
  evaluation.

  Can include an object with evaluator functions such as those generated
  by
  [`bru_eval_fun()`](https://inlabru-org.github.io/inlabru/reference/bru_eval_fun.md),
  names `<label>_eval`, one for each model component, `<label>`.

- pronouns, objects:

  If `NULL` (default), all named elements of `data` are used as
  pronouns/objects. If character vector(s), only the matching named
  elements of `data` will be available as pronouns/objects.

## Value

A data mask environment for use with
[`rlang::eval_tidy()`](https://rlang.r-lib.org/reference/eval_tidy.html).

## Examples

``` r
m <- bru_data_mask(
  data = list(
    data = data.frame(z = 11:14),
    extra_data = list(something = 1:2)
  )
)
rlang::eval_tidy(
  rlang::quo(z + rep(.extra_data$something, 2)),
  data = m
)
#> [1] 12 14 14 16

mask <- bru_data_mask(list(data = list(x = 1:4)))
rlang::eval_tidy(rlang::quo(x), data = mask)
#> [1] 1 2 3 4
rlang::eval_tidy(rlang::quo(.data$x), data = mask)
#> [1] 1 2 3 4
rlang::eval_tidy(rlang::quo(.data.$x), data = mask)
#> [1] 1 2 3 4

# Using functions that need to call eval_tidy using the data mask:
fun_factory <- function() {
  .get_mask <- function(frame = parent.frame(2L)) {
    rlang::eval_tidy(rlang::parse_expr(".mask."), env = frame)
  }
  .eval_tidy <- function(expr, frame = parent.frame(2L)) {
    mask <- .get_mask(frame)
    rlang::eval_tidy(expr, data = mask, env = frame)
  }
  .addition <- 5L
  fun <- function(nm) {
    val <- .eval_tidy(rlang::expr(.data[[nm]]))
    val + .addition
  }
  list(fun = fun)
}
mask <- bru_data_mask(list(data = list(x = 1:4), fun = fun_factory()))
rlang::eval_tidy(rlang::quo(fun("x")), data = mask)
#> [1] 6 7 8 9
rlang::eval_tidy(rlang::quo(.fun$fun("x")), data = mask)
#> [1] 6 7 8 9
rlang::eval_tidy(rlang::quo(.fun.$fun("x")), data = mask)
#> [1] 6 7 8 9
```
