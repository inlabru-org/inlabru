# Evaluate expressions in data contexts

Evaluate an expression in a series of data contexts, also making the
objects directly available as names surrounded by ".", stopping when the
expression evaluation completes with no error, as well as `tidy`
evaluation pronouns, e.g. `.data`.

This is an internal inlabru method, not intended for general use.

## Usage

``` r
bru_eval_in_data_context(
  input,
  data = NULL,
  default = NULL,
  .envir = parent.frame()
)
```

## Arguments

- input:

  An expression to be evaluated

- data:

  list of data objects in priority order. Named elements will be
  available as whole objects `.name.` as well as pronouns `.name` in the
  evaluation.

- default:

  Value used if the expression is evaluated as NULL. Default NULL

- .envir:

  The evaluation environment

## Value

The result of expression evaluation

## Examples

``` r
# The A values come from the 'data' element, and the B values come from
# the 'response_data' element, as that is listed first.
bru_eval_in_data_context(
  list(A = .data$x, B = x),
  list(
    response_data = tibble::tibble(x = 1:5),
    data = tibble::tibble(x = 1:10)
  )
)
#> $A
#>  [1]  1  2  3  4  5  6  7  8  9 10
#> 
#> $B
#> [1] 1 2 3 4 5
#> 
# Both A and B come from the 'data' element, as 'x' is found there,
# and the data objects are listed in order of precedence.
bru_eval_in_data_context(
  list(A = .data$x, B = x),
  list(
    data = tibble::tibble(x = 1:10),
    response_data = tibble::tibble(x = 1:5)
  )
)
#> $A
#>  [1]  1  2  3  4  5  6  7  8  9 10
#> 
#> $B
#>  [1]  1  2  3  4  5  6  7  8  9 10
#> 
```
