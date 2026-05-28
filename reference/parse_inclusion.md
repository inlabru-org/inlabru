# Parse inclusion of component labels in a predictor expression

Parse inclusion of component labels in a predictor expression

## Usage

``` r
parse_inclusion(thenames, include = NULL, exclude = NULL)
```

## Arguments

- thenames:

  Set of labels to restrict

- include:

  Character vector of component labels that are needed by the predictor
  expression; Default: NULL (include all components that are not
  explicitly excluded)

- exclude:

  Character vector of component labels that are not used by the
  predictor expression. The exclusion list is applied to the list as
  determined by the `include` parameter; Default: NULL (do not remove
  any components from the inclusion list)
