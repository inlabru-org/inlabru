# Transformation tools

Tools for transforming between N(0,1) variables and other distributions
in predictor expressions

## Usage

``` r
bru_forward_transformation(qfun, x, ..., tail.split. = 0)

bru_inverse_transformation(pfun, x, ..., tail.split. = NULL)
```

## Arguments

- qfun:

  A quantile function object, such as `qexp`

- x:

  Values to be transformed

- ...:

  Distribution parameters passed on to the `qfun` and `pfun` functions

- tail.split.:

  For x-values larger than `tail.split.`, upper quantile calculations
  are used internally, and for smaller values lower quantile
  calculations are used. This can avoid lack of accuracy in the
  distribution tails. If `NULL`, forward calculations split at 0, and
  inverse calculations use lower tails only, potentially losing accuracy
  in the upper tails.

- pfun:

  A CDF function object, such as `pexp`

## Value

- For `bru_forward_transformation`, a numeric vector

&nbsp;

- For `bru_inverse_transformation`, a numeric vector

## Examples

``` r
u <- rnorm(5, 0, 1)
y <- bru_forward_transformation(qexp, u, rate = 2)
v <- bru_inverse_transformation(pexp, y, rate = 2)
rbind(u, y, v)
#>          [,1]      [,2]     [,3]      [,4]     [,5]
#> u -1.78911563 0.2846580 1.194728 0.4832875 1.115947
#> y  0.01874611 0.4734355 1.076666 0.5784718 1.011635
#> v -1.78911563 0.2846580 1.194728 0.4832875 1.115947
```
