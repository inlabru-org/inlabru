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
#> 1  0.9584175 0.9437467 1.988361      0.03244409    0.04111276  5.589055
#> 2  1.0605426 1.0478110 1.877061      0.03587863    0.04338551  4.855778
#> 3  1.0497680 1.1237562 2.337244      0.03912182    0.05669211  8.178308
#> 4  1.0283020 0.9890142 1.622844      0.03344602    0.03432092  2.814957
#> 5  1.0017529 1.0169678 2.183948      0.03525910    0.04901147  7.288539
#> 6  0.9512145 0.9192355 1.595441      0.03105789    0.03145069  2.680380
#> 7  0.9963307 0.9981031 2.022379      0.03442999    0.04533432  6.250077
#> 8  1.0123914 0.9654978 1.701680      0.03281316    0.03607263  3.581580
#> 9  0.9760245 0.9621343 1.673310      0.03258573    0.03415843  3.039783
#> 10 1.0042888 1.0577224 2.308991      0.03676171    0.05239242  7.812165
```
