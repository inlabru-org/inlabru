# Interface between `bru_input` and `bru_mapper`

Associate
[bru_input](https://inlabru-org.github.io/inlabru/reference/bru_input.md)
objects with
[bru_mapper](https://inlabru-org.github.io/inlabru/reference/bru_mapper.md)
objects.

## Usage

``` r
ibm_input_set(mapper, input)

ibm_input_new(mapper, ...)

ibm_input_available(mapper)

ibm_input_get(mapper)

bm_autodetect()
```

## Arguments

- mapper:

  A `bru_mapper` object.

- input:

  A `bru_input` object, or `NULL` to remove any existing input.
  Alternatively, an existing `bru_mapper` object with or without
  associated input can be provided, in which case the input from that
  mapper is copied.

- ...:

  Passed on to
  [`new_bru_input()`](https://inlabru-org.github.io/inlabru/reference/bru_input.md).

## Functions

- `ibm_input_set()`: Add an existing `bru_input` to a `bru_mapper`.

- `ibm_input_new()`: Create and add a `bru_input` to a `bru_mapper`.

- `ibm_input_available()`: Check if a `bru_input` is associated with a
  `bru_mapper`.

- `ibm_input_get()`: Get the `bru_input` associated with a `bru_mapper`.

- `bm_autodetect()`: Create a `bru_mapper` placeholder object of class
  `bm_autodetect`. The main purpose of this class is to attach
  [bru_input](https://inlabru-org.github.io/inlabru/reference/bru_input.md)
  information to it, which is later used to determine a suitable 'real'
  mapper type.

## See also

[`bru_input()`](https://inlabru-org.github.io/inlabru/reference/bru_input.md)

## Examples

``` r
(m <- bm_autodetect())
#> autodetect
ibm_input_available(m)
#> [1] FALSE
(m <- ibm_input_new(m, cos(x)))
#> autodetect(cos(x))
ibm_input_available(m)
#> [1] TRUE
ibm_input_set(bm_linear(), m)
#> linear(cos(x))
```
