# Backwards compatibility to handle mexpand for INLA \<= 24.06.02

Expand observation vectors/matrices in stacks into to a multicolumn
matrix for multiple likelihoods

## Usage

``` r
bru_inla.stack.mexpand(
  ...,
  old.names = "BRU.response",
  new.name = "BRU.response"
)
```

## Arguments

- ...:

  List of stacks that contain vector observations (existing
  multilikelihood observation matrices are also permitted)

- old.names:

  A vector of strings with the names of the observation vector/matrix
  for each stack. If a single string, this is assumed for all the
  stacks. (default "BRU.response")

- new.name:

  The name to be used for the expanded observation matrix, possibly the
  same as an old name. (default "BRU.response")

## Value

A list of modified stacks with multicolumn observations

## Author

Fabian E. Bachl <f.e.bachl@bath.ac.uk> and Finn Lindgren
<finn.lindgren@gmail.com>
