# Methods for mapper extraction

Extract a mapper from another object

## Usage

``` r
as_bru_mapper(x)

# S3 method for class 'bru_mapper'
as_bru_mapper(x)

# S3 method for class 'bru_comp'
as_bru_mapper(x)

# S3 method for class 'bru_subcomp'
as_bru_mapper(x)
```

## Arguments

- x:

  Object to convert/extract

## Value

A `bru_mapper` object

## Examples

``` r
# Extract a mapper from a `bru_subcomp` object
as_bru_mapper(bru_comp("x", x, mapper = bm_index(4))$main)
#> index(x)
```
