# Extract timing information from fitted [bru](https://inlabru-org.github.io/inlabru/reference/bru.md) object

Extracts a data.frame or tibble with information about the `Time` (CPU),
`System`, and `Elapsed` time for each step of a
[`bru()`](https://inlabru-org.github.io/inlabru/reference/bru.md) run.

## Usage

``` r
bru_timings(object, ...)

# S3 method for class 'bru'
bru_timings(object, ...)
```

## Arguments

- object:

  A fitted `bru` object

- ...:

  unused

## Value

A `data.frame` or `tibble` with columns `Task`, `Iteration`, `Time`,
`System`, and `Elapsed`.
