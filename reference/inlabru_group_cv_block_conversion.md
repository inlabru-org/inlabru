# Convert node-block indices into groups of node indices

Converts a vector of block indices into a list of node index vectors.

## Usage

``` r
inlabru_group_cv_block_conversion(block, max_block, per_node)
```

## Arguments

- block:

  An integer vector of 1-based block indices, one for each node.

- max_block:

  Precomputed maximal value of `block`. Must be at least as large as
  `max(block)`. Block indices larger than `max_block` are ignored.

- per_node:

  Logical; if `FALSE`, the function returns a list with one element for
  each block. If `TRUE`, the function returns a list with one element
  for each node.

## Value

If `per_node` is `false`, a list of integer vectors, where each vector
contains the 1-based node indices of nodes belonging to the
corresponding block. If `per_node` is `true`, a list of integer vectors,
where each vector contains the 1-based node indices of nodes belonging
to the same block as the corresponding node; `b[i] <- a[block[j]]`
