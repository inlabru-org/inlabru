# Summarise and annotate data

Summarise and annotate data

## Usage

``` r
bru_summarise(
  data,
  probs = c(0.025, 0.5, 0.975),
  x = NULL,
  cbind.only = FALSE,
  max_moment = 2
)
```

## Arguments

- data:

  A list of samples, each either numeric or a `data.frame`

- probs:

  A numeric vector of probabilities with values in `[0, 1]`, passed to
  [`stats::quantile`](https://rdrr.io/r/stats/quantile.html)

- x:

  A `data.frame` of data columns that should be added to the summary
  data frame

- cbind.only:

  If TRUE, only `cbind` the samples and return a matrix where each
  column is a sample

- max_moment:

  integer, at least 2. Determines the largest moment order information
  to include in the output. If `max_moment > 2`, includes "skew"
  (skewness, `E[(x-m)^3/s^3]`), and if `max_moment > 3`, includes
  "ekurtosis" (excess kurtosis, `E[(x-m)^4/s^4] - 3`). Default 2. Note
  that the Monte Carlo variability of the `ekurtois` estimate may be
  large.

## Value

A `data.frame` or `Spatial[Points/Pixels]DataFrame` with summary
statistics, "mean", "sd", `paste0("q", probs)`, "mean.mc_std_err",
"sd.mc_std_err"

## Examples

``` r
bru_summarise(matrix(rexp(10000), 10, 1000), max_moment = 4, probs = NULL)
#>         mean        sd     skew mean.mc_std_err sd.mc_std_err ekurtosis
#> 1  0.9558905 0.9389345 2.009806      0.03230460    0.04131330  5.742065
#> 2  1.0611409 1.0489573 1.873341      0.03591258    0.04334911  4.829320
#> 3  1.0536426 1.1239889 2.331875      0.03912434    0.05661570  8.146682
#> 4  1.0299244 0.9892172 1.623983      0.03345313    0.03433187  2.816050
#> 5  1.0037306 1.0177351 2.168777      0.03527648    0.04890251  7.233336
#> 6  0.9526178 0.9267548 1.623076      0.03133645    0.03209540  2.795507
#> 7  0.9973735 0.9967847 2.027467      0.03439063    0.04537118  6.285380
#> 8  1.0132720 0.9639491 1.702073      0.03276525    0.03608960  3.604808
#> 9  0.9733775 0.9625180 1.674963      0.03259972    0.03418782  3.044439
#> 10 1.0068927 1.0551332 2.316476      0.03668425    0.05246231  7.886724
```
