# Extract basic variable names from expression

Extracts the variable names from an R expression by pre- and
post-processing around
[`all.vars()`](https://rdrr.io/r/base/allnames.html). First replaces `$`
with `[[` indexing, so that internal column/variable names are ignored,
then calls [`all.vars()`](https://rdrr.io/r/base/allnames.html).

## Usage

``` r
bru_used_vars(x, functions = FALSE)

# S3 method for class 'character'
bru_used_vars(x, functions = FALSE)

# S3 method for class 'expression'
bru_used_vars(x, functions = FALSE)

# S3 method for class 'quosure'
bru_used_vars(x, functions = FALSE)

# S3 method for class 'formula'
bru_used_vars(x, functions = FALSE)
```

## Arguments

- x:

  A `formula`, `expression`, or `character`

- functions:

  logical; if TRUE, include function names

## Value

If successful, a character vector, otherwise `NULL`

## Methods (by class)

- `bru_used_vars(formula)`: Only the right-hand side is used.

## See also

Other bru_used:
[`bru_used()`](https://inlabru-org.github.io/inlabru/reference/bru_used.md),
[`bru_used_update()`](https://inlabru-org.github.io/inlabru/reference/bru_used_update.md),
[`new_bru_used()`](https://inlabru-org.github.io/inlabru/reference/new_bru_used.md)

## Examples

``` r
bru_used_vars(~.)
#> NULL
bru_used_vars(~ a + b + c_latent + d_eval())
#> [1] "a"        "b"        "c_latent"
bru_used_vars(expression(a + b + c_latent + d_eval()))
#> [1] "a"        "b"        "c_latent"

bru_used_vars(~., functions = TRUE)
#> NULL
bru_used_vars(~ a + b + c_latent + d_eval(), functions = TRUE)
#> [1] "+"        "a"        "b"        "c_latent" "d_eval"  
bru_used_vars(expression(a + b + c_latent + d_eval()), functions = TRUE)
#> [1] "expression" "+"          "a"          "b"          "c_latent"  
#> [6] "d_eval"    

bru_used_vars(a ~ b)
#> [1] "b"
bru_used_vars(expression(a ~ b))
#> [1] "a" "b"
```
