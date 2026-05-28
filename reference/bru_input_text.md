# Extract input definitions as text

Extract input definitions from the
[bru_input](https://inlabru-org.github.io/inlabru/reference/bru_input.md)
information stored in bru model components. Primarily intended for
diagnostic use, as the output does not currently contain all
information, such as the layer and selector information for spatial
covariates.

## Usage

``` r
bru_input_text(...)

# S3 method for class 'bru_input'
bru_input_text(x, .envir = parent.frame(), ...)

# S3 method for class 'bru_comp'
bru_input_text(x, ..., label = x$label)

# S3 method for class 'bru_comp_list'
bru_input_text(x, ...)

# S3 method for class 'bru_mapper'
bru_input_text(x, ..., label = "<unknown>")

# S3 method for class 'bm_pipe'
bru_input_text(x, ..., label = "<unknown>")

# S3 method for class 'bm_multi'
bru_input_text(x, ..., label = "<unknown>")

# S3 method for class 'bm_collect'
bru_input_text(x, ..., label = "<unknown>")

# S3 method for class 'bm_repeat'
bru_input_text(x, ..., label = "<unknown>")

# S3 method for class 'bm_sum'
bru_input_text(x, ..., label = "<unknown>")
```

## Arguments

- ...:

  Passed on to sub-methods.

- x:

  A `bru_input` object, or other object for recursive evaluation.

- .envir:

  environment in which to evaluate the input expression. Default is
  [`parent.frame()`](https://rdrr.io/r/base/sys.parent.html)

- label:

  character; optional label used to identify the object in informational
  messages.

## Value

- `bru_input(bru_comp)`: A character vector of mapper input values.

&nbsp;

- `bru_input_text(bru_comp_list)`: A list of mapper input definition
  text strings, with one entry for each component.

## Methods (by class)

- `bru_input_text(bru_input)`: Extract the input definitions from a
  `bru_input` object, as a named character vector

- `bru_input_text(bru_comp)`: Extract the input definition as a named
  character vector.

- `bru_input_text(bru_mapper)`: Extract the input definitions associated
  with a `bru_mapper`.

## See also

[`bru_input()`](https://inlabru-org.github.io/inlabru/reference/bru_input.md)

[ibm_input](https://inlabru-org.github.io/inlabru/reference/ibm_input.md)

## Examples

``` r
(inp <- new_bru_input(x, "LABEL"))
#> LABEL = x
bru_input_text(inp)
#> [1] "x"
if (bru_safe_inla()) {
  bru_input_text(bru_comp("x", cos(y)))
}
#> $core
#> $core$main
#> [1] "cos(y)"
#> 
#> 
```
