# Re-indexing matrix for repeated mappers

Helper methods for regular and interleaved state vector indexing, meant
for use in
[bm_repeat](https://inlabru-org.github.io/inlabru/reference/bm_repeat.md)
mappers.

## Usage

``` r
bm_repeat_indexing(n_map, n_rep, interleaved = FALSE)

bm_repeat_indexing_matrix(n_map, n_rep)
```

## Arguments

- n_map:

  The block mapper size

- n_rep:

  The number of times the block is repeated

- interleaved:

  logical; if `TRUE`, the state vector indexing is interleaved. Default
  is `FALSE`, for blockwise indexing.

## Value

`bm_repeat_indexing`: Returns a list with a vector of offsets,
`offsets`, and a vector of relative index values, `index`; block `k`
indexing is given by `offsets[k] + index`.

`bm_repeat_indexing_matrix`: A `sparseMatrix` object.

## Functions

- `bm_repeat_indexing()`: Construct the offsets and within-block index
  vectors for a repeated mapper, with or without interleaving.

- `bm_repeat_indexing_matrix()`: Creates a sparse matrix `A` such that
  `z <- A %*% x` constructs a blockwise version `z = c(x1, x2, ..., xn)`
  of an interleaved state vector
  `x = c(x1[1], x2[1], ..., x1[2], x2[2], ...)`. Each block is of size
  `n_map`, and there are `n_rep` blocks. The reverse operation, taking a
  blockwise `z = c(x1, x2, ..., xn)` to an interleaved vector is
  `x <- Matrix::t(A) %*% z`.

## Examples

``` r
(idx <- bm_repeat_indexing(3, 2, FALSE))
#> $offsets
#> [1] 0 3
#> 
#> $index
#> [1] 1 2 3
#> 
(idx <- bm_repeat_indexing(3, 2, TRUE))
#> $offsets
#> [1] 0 1
#> 
#> $index
#> [1] 1 3 5
#> 
(A <- bm_repeat_indexing_matrix(3, 2))
#> 6 x 6 sparse Matrix of class "dgCMatrix"
#>                 
#> [1,] 1 . . . . .
#> [2,] . . 1 . . .
#> [3,] . . . . 1 .
#> [4,] . 1 . . . .
#> [5,] . . . 1 . .
#> [6,] . . . . . 1
(x_interleaved <- 1:6)
#> [1] 1 2 3 4 5 6
(x_blockwise <- as.vector(A %*% x_interleaved))
#> [1] 1 3 5 2 4 6
(x_recovered <- as.vector(Matrix::t(A) %*% x_blockwise))
#> [1] 1 2 3 4 5 6
```
