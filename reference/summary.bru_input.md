# Summarise component inputs

Summarise component inputs

## Usage

``` r
# S3 method for class 'bru_input'
format(x, verbose = TRUE, ..., label.override = NULL, type = NULL)

# S3 method for class 'bru_input'
summary(object, verbose = TRUE, ..., label.override = NULL)

# S3 method for class 'bru_input'
print(x, verbose = TRUE, ..., label.override = NULL)
```

## Arguments

- x:

  Object to be printed

- verbose:

  logical; If `TRUE`, includes more details of the component
  definitions. When `FALSE`, only show basic component definition
  information. Default `TRUE`.

- ...:

  Passed on to other summary methods.

- label.override:

  character; If not `NULL`, use this label instead of the object's
  label.

- type:

  character; if non-NULL, added to the output'; `label = type(input)`.

- object:

  Object to be summarised.

## See also

[`bru_input()`](https://inlabru-org.github.io/inlabru/reference/bru_input.md),
[`bru_comp()`](https://inlabru-org.github.io/inlabru/reference/bru_comp.md)

## Author

Fabian E. Bachl <bachlfab@gmail.com>

Finn Lindgren <finn.lindgren@gmail.com>
