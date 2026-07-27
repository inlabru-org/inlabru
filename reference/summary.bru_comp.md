# Summarise components

Summarise components

## Usage

``` r
# S3 method for class 'summary_bru_input'
print(x, ...)

# S3 method for class 'bru_comp'
summary(object, ..., depth = Inf, verbose = TRUE)

# S3 method for class 'bru_comp_list'
summary(object, verbose = TRUE, ...)

# S3 method for class 'bru_comp'
print(x, ...)

# S3 method for class 'bru_comp_list'
print(x, ...)

# S3 method for class 'bru_subcomp'
format(x, verbose = TRUE, ..., label.override = NULL)

# S3 method for class 'bru_subcomp'
summary(object, verbose = TRUE, ..., label.override = NULL)

# S3 method for class 'bru_subcomp'
print(x, verbose = TRUE, ..., label.override = NULL)

# S3 method for class 'summary_bru_comp'
print(x, ...)

# S3 method for class 'summary_bru_comp_list'
print(x, ...)

# S3 method for class 'summary_bru_subcomp'
print(x, ...)
```

## Arguments

- x:

  A summary object to be printed.

- ...:

  Passed on to other summary methods.

- object:

  Object to be summarised.

- depth:

  The depth of which to expand the component mapper. Default `Inf`, to
  traverse the entire mapper tree.

- verbose:

  logical; If `TRUE`, includes more details of the component
  definitions. When `FALSE`, only show basic component definition
  information. Default `TRUE`.

- label.override:

  character; If not `NULL`, use this label instead of the object's
  label.

## See also

[`bru_comp()`](https://inlabru-org.github.io/inlabru/reference/bru_comp.md),
[`summary.bru_input()`](https://inlabru-org.github.io/inlabru/reference/summary.bru_input.md)

## Author

Fabian E. Bachl <bachlfab@gmail.com>

Finn Lindgren <finn.lindgren@gmail.com>
