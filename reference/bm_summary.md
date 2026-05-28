# mapper object summaries

mapper object summaries

## Usage

``` r
# S3 method for class 'bru_mapper'
format(x, ..., prefix = "", initial = prefix, depth = 1)

# S3 method for class 'bm_list'
format(
  x,
  ...,
  prefix = "",
  initial = prefix,
  depth = 1,
  collapse = ", ",
  labels = TRUE
)

# S3 method for class 'bru_mapper'
summary(object, ..., prefix = "", initial = prefix, depth = 1)

# S3 method for class 'bm_multi'
format(x, ..., prefix = "", initial = prefix, depth = 1)

# S3 method for class 'bm_pipe'
format(x, ..., prefix = "", initial = prefix, depth = 1)

# S3 method for class 'bm_collect'
format(x, ..., prefix = "", initial = prefix, depth = 1)

# S3 method for class 'bm_sum'
format(x, ..., prefix = "", initial = prefix, depth = 1)

# S3 method for class 'bm_repeat'
format(x, ..., prefix = "", initial = prefix, depth = 1)

# S3 method for class 'bm_reparam'
format(x, ..., prefix = "", initial = prefix, depth = 1)

# S3 method for class 'summary_bru_mapper'
print(x, ..., sep = "\n")

# S3 method for class 'bru_mapper'
print(x, ..., sep = "\n", prefix = "", initial = prefix, depth = 1)

# S3 method for class 'bm_list'
print(
  x,
  ...,
  sep = "\n",
  prefix = "",
  initial = prefix,
  depth = 1,
  labels = TRUE,
  collapse = ", "
)
```

## Arguments

- x:

  Object to format/print

- ...:

  Unused arguments

- prefix:

  character prefix for each line. Default `""`.

- initial:

  character prefix for the first line. Default `initial=prefix`.

- depth:

  The recursion depth for multi/collection/pipe mappers. Default 1, to
  only show the collection, and not the contents of the sub-mappers.

- collapse:

  character or NULL, as in
  [`base::paste()`](https://rdrr.io/r/base/paste.html).

- labels:

  logical; if TRUE, include mapper names or numerical indices. Default
  `TRUE`

- object:

  Object to summarise

- sep:

  character; separator for printing the summary.

## Examples

``` r
mapper <-
  bm_pipe(
    list(
      bm_multi(list(
        A = bm_index(2),
        B = bm_index(3)
      )),
      bm_index(2)
    )
  )
summary(mapper, depth = 2)
#> pipe = multi(A = index, B = index) -> index
mapper <-
  bm_repeat(
    bm_multi(
      list(
        A = bm_index(2),
        B = bm_index(3)
      )
    ),
    3
  )
summary(mapper)
#> repeat(3 x multi(A = index, B = index))
summary(mapper, depth = 0)
#> repeat(3 x multi)
mapper <-
  bm_reparam(
    bm_multi(
      list(
        A = bm_index(2),
        B = bm_index(3)
      )
    ),
    matrix(1:36, nrow = 6)
  )
summary(mapper)
#> reparam(multi(A = index, B = index))
summary(mapper, depth = 0)
#> reparam(multi)
```
