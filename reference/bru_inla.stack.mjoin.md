# Join stacks intended to be run with different likelihoods

Helper functions for multi-likelihood models

## Usage

``` r
bru_inla.stack.mjoin(
  ...,
  compress = TRUE,
  remove.unused = TRUE,
  old.names = "BRU.response",
  new.name = "BRU.response"
)
```

## Arguments

- ...:

  List of stacks that contain vector observations (existing
  multi-likelihood observation matrices are also permitted)

- compress:

  If `TRUE`, compress the model by removing duplicated rows of effects,
  replacing the corresponding A-matrix columns with a single column
  containing the sum.

- remove.unused:

  If `TRUE`, compress the model by removing rows of effects
  corresponding to all-zero columns in the A matrix (and removing those
  columns).

- old.names:

  A vector of strings with the names of the observation vector/matrix
  for each stack. If a single string, this is assumed for all the
  stacks. (default "BRU.response")

- new.name:

  The name to be used for the expanded observation matrix, possibly the
  same as an old name. (default "BRU.response")

## Value

A single stack with a multi-likelihood observation matrix
