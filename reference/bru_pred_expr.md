# Predictor expression handler object

Create and access predictor expression information, keeping track of
what latent variables are used, whether the predictor is
additive/linear/rowwise.

## Usage

``` r
new_bru_pred_expr(
  x,
  ...,
  used = NULL,
  is_rowwise = NULL,
  .envir = parent.frame()
)

bru_pred_expr(x, ...)

# S3 method for class 'bru_obs'
bru_pred_expr(x, ...)

# S3 method for class 'bru_obs_list'
bru_pred_expr(x, ...)

# S3 method for class 'bru_info'
bru_pred_expr(x, ...)

# S3 method for class 'bru'
bru_pred_expr(x, ...)

# S3 method for class 'bru_pred_expr'
bru_pred_expr(x, ..., format = "object")

# S3 method for class 'bru_pred_expr'
format(x, ...)

# S3 method for class 'bru_pred_expr'
print(x, ...)
```

## Arguments

- x:

  For creation, a formula or character string giving the predictor
  expression. If the character string is `"."`, it is taken to mean an
  additive expression including all latent variables. For access, an
  object containing a `bru_pred_expr` object.

- ...:

  Passed on to submethods

- used:

  An optional `bru_used` object, overriding the automated detection.
  Default: `NULL`

- is_rowwise:

  logical; whether the predictor can be assumed to use the input data
  rowwise, so that blockwise evaluation is possible. If `NULL`, it needs
  to be set before use. Default: `NULL`

- .envir:

  The environment in which to evaluate the expression.

- format:

  Character; one of

  - `"object"` (the default),

  - `"text"` (plain text),

  - `"quo"` (an `rlang` quosure),

  - `"expr"` (an `rlang` expression),

  - `"formula"` (the original formula input if available, otherwise
    constructed from the expression),

  - `"formula_text"` (a text version of the "formula"),

  - `"text_raw"` (plain text, without substituting `BRU_EXPRESSION`), or

  - `"resp_text"` (plain text of the response side of the formula).

## Methods (by class)

- `bru_pred_expr(bru_obs)`: Accessor for the `bru_pred_expr` object
  stored inside a `bru_obs` object.

- `bru_pred_expr(bru_obs_list)`: Accessor for the `bru_pred_expr`
  objects stored inside a `bru_obs_list` object.

- `bru_pred_expr(bru_info)`: Accessor for the `bru_pred_expr` object
  stored inside a `bru_info` object.

- `bru_pred_expr(bru)`: Accessor for the `bru_pred_expr` object stored
  inside a `bru` object.

- `bru_pred_expr(bru_pred_expr)`: Access a `bru_pred_expr` object or
  convert it to text or expression object.

## Functions

- `new_bru_pred_expr()`: Create a `bru_pred_expr` object from a formula
  or character string.

- `bru_pred_expr()`: Accessor generic for `bru_pred_expr` objects,
  including ones stored inside other objects.

## Examples

``` r
(new_bru_pred_expr(~ x + z + Intercept))
#> Predictor: ~ x + z + Intercept
#> Additive/Linear/Rowwise: TRUE/TRUE/FALSE
#> Used components: effect[x, z, Intercept], latent[]
```
