# Extract basic variable names from expression

Extracts the variable names and function names from an R expression by
traversing the expression structure. Internal helper function for
[`new_bru_used()`](https://inlabru-org.github.io/inlabru/reference/new_bru_used.md)
and
[`bru_used()`](https://inlabru-org.github.io/inlabru/reference/bru_used.md).

## Usage

``` r
bru_used_vars(x, result = new_bru_used_vars())

new_bru_used_vars(
  x = list(vars = character(0), funs = character(0), objects = list())
)

# Default S3 method
bru_used_vars(x, result = new_bru_used_vars())

# S3 method for class '`<-`'
bru_used_vars(x, result = new_bru_used_vars())

# S3 method for class 'call'
bru_used_vars(x, result = new_bru_used_vars())

# S3 method for class 'expression'
bru_used_vars(x, result = new_bru_used_vars())

# S3 method for class 'quosure'
bru_used_vars(x, result = new_bru_used_vars())

# S3 method for class 'formula'
bru_used_vars(x, result = new_bru_used_vars())

# S3 method for class 'bru_used_vars'
format(x, ...)

# S3 method for class 'bru_used_vars'
print(x, ...)
```

## Arguments

- x:

  A `formula`, `expression`, or other supported class. For the `format`
  and `print` methods, a `bru_used_vars` object.

- result:

  A `bru_used_vars` object; a list with elements `vars`, `funs`, and
  `objects`, by default provided by `new_bru_used_vars()`

## Value

A `bru_used_vars` object with elements

- vars:

  character; names of directly accessed variables.

- funs:

  character; names of functions called.

- objects:

  named list; one character vector per container objects with variables
  names accessed via `$`, `[[`, or `[`. If the access is ambiguous, the
  container object name is stored in `vars`.

## Functions

- `new_bru_used_vars()`: Create a `bru_used_vars` object.

## See also

Other bru_used:
[`bru_used()`](https://inlabru-org.github.io/inlabru/reference/bru_used.md),
[`bru_used_update()`](https://inlabru-org.github.io/inlabru/reference/bru_used_update.md),
[`new_bru_used()`](https://inlabru-org.github.io/inlabru/reference/new_bru_used.md)

## Examples

``` r
bru_used_vars(~.)
#> vars: {.}, funs: {}, objects[]
bru_used_vars(~ a + b + c_latent + d_eval())
#> vars: {a, b, c_latent}, funs: {+, d_eval}, objects[]

# Ignores the LHS:
bru_used_vars(a ~ b)
#> vars: {b}, funs: {}, objects[]

# Detects variables accessed via pronouns and objects,
# as well as function calls:
bru_used_vars(~ cos(x$z) + y_eval() + .latent$"q")
#> vars: {}, funs: {+, cos, y_eval}, objects[x: {z}, .latent: {q}]
```
