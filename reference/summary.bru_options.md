# Print inlabru options

Print inlabru options

## Usage

``` r
# S3 method for class 'bru_options'
summary(
  object,
  legend = TRUE,
  include_global = TRUE,
  include_default = TRUE,
  ...
)

# S3 method for class 'summary_bru_options'
print(x, ...)
```

## Arguments

- object:

  A
  [bru_options](https://inlabru-org.github.io/inlabru/reference/bru_options.md)
  object to be summarised

- legend:

  logical; If `TRUE`, include explanatory text, Default: `TRUE`

- include_global:

  logical; If `TRUE`, include global override options

- include_default:

  logical; If `TRUE`, include default options

- ...:

  Further parameters, currently ignored

- x:

  A `summary_bru_options` object to be printed

## Examples

``` r
if (interactive()) {
  options <- bru_options(verbose = TRUE)

  # Don't print options only set in default:
  print(options, include_default = FALSE)

  # Only include options set in the object:
  print(options, include_default = FALSE, include_global = FALSE)
}
```
