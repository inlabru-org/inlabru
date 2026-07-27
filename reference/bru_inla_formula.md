# Extract INLA formula

Extracts the INLA formula used in a `bru` model

## Usage

``` r
bru_inla_formula(x, ...)

# S3 method for class 'bru'
bru_inla_formula(x, ...)

# S3 method for class 'bru_model'
bru_inla_formula(x, ...)

# S3 method for class 'bru_info'
bru_inla_formula(x, ...)

# S3 method for class 'bru_comp_list'
bru_inla_formula(x, ...)
```

## Arguments

- x:

  An object containing information about an INLA formula

- ...:

  Additional arguments passed on to submethods

## Value

A [formula](https://rdrr.io/r/stats/formula.html) suited for use in
[`INLA::inla()`](https://rdrr.io/pkg/INLA/man/inla.html)
