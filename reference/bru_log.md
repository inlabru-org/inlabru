# Access methods for `bru_log` objects

Access method for `bru_log` objects. Note: Up to version `2.8.0`,
`bru_log()` was a deprecated alias for
[`bru_log_message()`](https://inlabru-org.github.io/inlabru/reference/bru_log_message.md).
When running on `2.8.0` or earlier, use `bru_log_get()` to access the
global log, and `cat(fit$bru_iinla$log, sep = "\n")` to print a stored
estimation object log. After version `2.8.0`, use `bru_log()` to access
the global log, and `bru_log(fit)` to access a stored estimation log.

## Usage

``` r
bru_log(x = NULL, verbosity = NULL)

# S3 method for class 'character'
bru_log(x, verbosity = NULL)

# S3 method for class 'bru_log'
bru_log(x, verbosity = NULL)

# S3 method for class 'iinla'
bru_log(x, verbosity = NULL)

# S3 method for class 'bru'
bru_log(x, verbosity = NULL)

# S3 method for class 'bru_log'
format(x, ..., timestamp = TRUE, verbosity = FALSE)

# S3 method for class 'bru_log'
print(x, ..., timestamp = TRUE, verbosity = FALSE)

# S3 method for class 'bru_log'
as.character(x, ...)

# S3 method for class 'bru_log'
x[i]

# S3 method for class 'bru_log'
c(...)

# S3 method for class 'bru_log'
length(x)
```

## Arguments

- x:

  An object that is, contains, or can be converted to, a `bru_log`
  object. If `NULL`, refers to the global `inlabru` log.

- verbosity:

  integer value for limiting the highest verbosity level being returned.

- ...:

  further arguments passed to or from other methods.

- timestamp:

  If `TRUE`, include the timestamp of each message. Default `TRUE`.

- i:

  indices specifying elements to extract. If `character`, denotes the
  sequence between bookmark `i` and the next bookmark (or the end of the
  log if `i` is the last bookmark)

## Value

`bru_log` A `bru_log` object, containing a character vector of log
messages, and potentially a vector of bookmarks.

## Methods (by generic)

- `format(bru_log)`: Format a `bru_log` object for printing. If
  `verbosity` is `TRUE`, include the verbosity level of each message.

- `print(bru_log)`: Print a `bru_log` object with `cat(x, sep = "\n")`.
  If `verbosity` is `TRUE`, include the verbosity level of each message.

- `as.character(bru_log)`: Convert `bru_log` object to a plain
  `character` vector

- `[`: Extract a subset of a `bru_log` object

- `c(bru_log)`: Concatenate several `bru_log` or `character` objects
  into a `bru_log` object.

- `length(bru_log)`: Obtain the number of log entries into a `bru_log`
  object.

## Functions

- `bru_log()`: Extract stored log messages. If non-`NULL`, the
  `verbosity` argument determines the maximum verbosity level of the
  messages to extract.

## See also

Other inlabru log methods:
[`bru_log_bookmark()`](https://inlabru-org.github.io/inlabru/reference/bru_log_bookmark.md),
[`bru_log_message()`](https://inlabru-org.github.io/inlabru/reference/bru_log_message.md),
[`bru_log_offset()`](https://inlabru-org.github.io/inlabru/reference/bru_log_offset.md),
[`bru_log_reset()`](https://inlabru-org.github.io/inlabru/reference/bru_log_reset.md),
[`new_bru_log()`](https://inlabru-org.github.io/inlabru/reference/new_bru_log.md)

## Examples

``` r
bru_log(verbosity = 2L)
#> 2026-07-27 22:31:52.447286: inlabru loaded
#> 2026-07-27 22:31:52.44774: Clear override options
#> 2026-07-27 22:32:05.535344: bru: Preprocessing
#> 2026-07-27 22:32:05.65072: iinla: Iteration 1 [max: 1]
#> 2026-07-27 22:32:06.623265: bru: Preprocessing
#> 2026-07-27 22:32:06.709365: iinla: Iteration 1 [max: 1]
#> 2026-07-27 22:32:07.205802: bru: Preprocessing
#> 2026-07-27 22:32:07.303859: iinla: Iteration 1 [max: 10]
#> 2026-07-27 22:32:07.765726: iinla: Step rescaling: 27.4% (norm0 = 184.1, norm1 = 225.5, norm01 = 347.8)
#> 2026-07-27 22:32:07.785129: iinla: Iteration 2 [max: 10]
#> 2026-07-27 22:32:08.543089: iinla: Step rescaling: 99.7% (norm0 = 224.6, norm1 = 10.68, norm01 = 225.4)
#> 2026-07-27 22:32:08.561785: iinla: Max deviation from previous: 46300% of SD, and line search is active
#> [stop if: < 10% and line search inactive]
#> 2026-07-27 22:32:08.563622: iinla: Iteration 3 [max: 10]
#> 2026-07-27 22:32:08.98853: iinla: Step rescaling: 102% (norm0 = 10.68, norm1 = 0.01176, norm01 = 10.68)
#> 2026-07-27 22:32:09.007997: iinla: Max deviation from previous: 496% of SD, and line search is active
#> [stop if: < 10% and line search inactive]
#> 2026-07-27 22:32:09.009904: iinla: Iteration 4 [max: 10]
#> 2026-07-27 22:32:09.475161: iinla: Max deviation from previous: 8.05% of SD, and line search is inactive
#> [stop if: < 10% and line search inactive]
#> 2026-07-27 22:32:09.476402: iinla: Convergence criterion met.
#>        Running final INLA integration step with known theta mode.
#> 2026-07-27 22:32:09.478259: iinla: Iteration 5 [max: 10]
#> 2026-07-27 22:32:11.836655: bru: Preprocessing
#> 2026-07-27 22:32:11.895514: iinla: Iteration 1 [max: 1]
#> 2026-07-27 22:32:12.455093: bru: Preprocessing
#> 2026-07-27 22:32:12.53919: iinla: Iteration 1 [max: 1]
#> 2026-07-27 22:32:14.953495: bru: Preprocessing
#> 2026-07-27 22:32:15.032713: iinla: Iteration 1 [max: 1]
#> 2026-07-27 22:32:23.209219: bru: Preprocessing
#> 2026-07-27 22:32:23.265291: iinla: Iteration 1 [max: 1]
#> 2026-07-27 22:32:25.921234: bru: Preprocessing
#> 2026-07-27 22:32:26.12987: iinla: Iteration 1 [max: 1]
#> 2026-07-27 22:33:09.761333: bru: Preprocessing
format(bru_log())
#> 2026-07-27 22:31:52.447286: inlabru loaded
#> 2026-07-27 22:31:52.44774: Clear override options
#> 2026-07-27 22:32:05.535344: bru: Preprocessing
#> 2026-07-27 22:32:05.557014: Evaluate component inputs for each observation model
#> 2026-07-27 22:32:05.559441: bru_input(bru_comp_list)
#> 2026-07-27 22:32:05.561167: bru_input.bru_comp(x)
#> 2026-07-27 22:32:05.562665: bru_input.bm_pipe(x)
#> 2026-07-27 22:32:05.564245: bru_input.bm_multi(x:core)
#> 2026-07-27 22:32:05.565727: bru_input.bru_mapper(x:core:main)
#> 2026-07-27 22:32:05.567357: bru_input(bru_input) for (x)
#> 2026-07-27 22:32:05.574951: bru_input.bru_comp(Intercept)
#> 2026-07-27 22:32:05.576444: bru_input.bm_pipe(Intercept)
#> 2026-07-27 22:32:05.578092: bru_input.bm_multi(Intercept:core)
#> 2026-07-27 22:32:05.579476: bru_input.bru_mapper(Intercept:core:main)
#> 2026-07-27 22:32:05.580757: bru_input(bru_input) for (Intercept)
#> 2026-07-27 22:32:05.596993: iinla: Start
#> 2026-07-27 22:32:05.59862: iinla: Evaluate component simplifications
#> 2026-07-27 22:32:05.599805: Simplify component mappers for each observation model
#> 2026-07-27 22:32:05.601259: Simplify component 'NA'
#> 2026-07-27 22:32:05.604469: Simplify component 'NA'
#> 2026-07-27 22:32:05.607868: iinla: Evaluate predictor linearisation
#> 2026-07-27 22:32:05.60991: Linearise predictor for each observation model
#> 2026-07-27 22:32:05.630034: iinla: Construct inla stack
#> 2026-07-27 22:32:05.649276: iinla: Model initialisation completed
#> 2026-07-27 22:32:05.65072: iinla: Iteration 1 [max: 1]
#> 2026-07-27 22:32:06.600036: iinla: Computation completed
#> 2026-07-27 22:32:06.623265: bru: Preprocessing
#> 2026-07-27 22:32:06.632179: Evaluate component inputs for each observation model
#> 2026-07-27 22:32:06.633593: bru_input(bru_comp_list)
#> 2026-07-27 22:32:06.634923: bru_input.bru_comp(x)
#> 2026-07-27 22:32:06.636266: bru_input.bm_pipe(x)
#> 2026-07-27 22:32:06.63777: bru_input.bm_multi(x:core)
#> 2026-07-27 22:32:06.639199: bru_input.bru_mapper(x:core:main)
#> 2026-07-27 22:32:06.640574: bru_input(bru_input) for (x)
#> 2026-07-27 22:32:06.648298: bru_input.bru_comp(Intercept)
#> 2026-07-27 22:32:06.649708: bru_input.bm_pipe(Intercept)
#> 2026-07-27 22:32:06.651115: bru_input.bm_multi(Intercept:core)
#> 2026-07-27 22:32:06.652484: bru_input.bru_mapper(Intercept:core:main)
#> 2026-07-27 22:32:06.653812: bru_input(bru_input) for (Intercept)
#> 2026-07-27 22:32:06.666989: iinla: Start
#> 2026-07-27 22:32:06.668622: iinla: Evaluate component simplifications
#> 2026-07-27 22:32:06.669775: Simplify component mappers for each observation model
#> 2026-07-27 22:32:06.671082: Simplify component 'NA'
#> 2026-07-27 22:32:06.674181: Simplify component 'NA'
#> 2026-07-27 22:32:06.678353: iinla: Evaluate predictor linearisation
#> 2026-07-27 22:32:06.681037: Linearise predictor for each observation model
#> 2026-07-27 22:32:06.690596: iinla: Construct inla stack
#> 2026-07-27 22:32:06.707399: iinla: Model initialisation completed
#> 2026-07-27 22:32:06.709365: iinla: Iteration 1 [max: 1]
#> 2026-07-27 22:32:07.178737: iinla: Computation completed
#> 2026-07-27 22:32:07.205802: bru: Preprocessing
#> 2026-07-27 22:32:07.216748: Evaluate component inputs for each observation model
#> 2026-07-27 22:32:07.218695: bru_input(bru_comp_list)
#> 2026-07-27 22:32:07.220516: bru_input.bru_comp(z)
#> 2026-07-27 22:32:07.222291: bru_input.bm_pipe(z)
#> 2026-07-27 22:32:07.224382: bru_input.bm_multi(z:core)
#> 2026-07-27 22:32:07.226321: bru_input.bru_mapper(z:core:main)
#> 2026-07-27 22:32:07.228101: bru_input(bru_input) for (z)
#> 2026-07-27 22:32:07.235719: bru_input.bru_comp(Intercept)
#> 2026-07-27 22:32:07.237532: bru_input.bm_pipe(Intercept)
#> 2026-07-27 22:32:07.239467: bru_input.bm_multi(Intercept:core)
#> 2026-07-27 22:32:07.241391: bru_input.bru_mapper(Intercept:core:main)
#> 2026-07-27 22:32:07.243237: bru_input(bru_input) for (Intercept)
#> 2026-07-27 22:32:07.258339: iinla: Start
#> 2026-07-27 22:32:07.260482: iinla: Evaluate component simplifications
#> 2026-07-27 22:32:07.262089: Simplify component mappers for each observation model
#> 2026-07-27 22:32:07.263993: Simplify component 'NA'
#> 2026-07-27 22:32:07.26787: Simplify component 'NA'
#> 2026-07-27 22:32:07.272138: iinla: Evaluate predictor linearisation
#> 2026-07-27 22:32:07.274735: Linearise predictor for each observation model
#> 2026-07-27 22:32:07.284394: iinla: Construct inla stack
#> 2026-07-27 22:32:07.301619: iinla: Model initialisation completed
#> 2026-07-27 22:32:07.303859: iinla: Iteration 1 [max: 10]
#> 2026-07-27 22:32:07.757508: iinla: Step rescaling: 61.8%, Contract (norm0 = 1969, norm1 = 1807, norm01 = 347.8)
#> 2026-07-27 22:32:07.759907: iinla: Step rescaling: 38.2%, Contract (norm0 = 398.9, norm1 = 295.1, norm01 = 347.8)
#> 2026-07-27 22:32:07.762973: iinla: Step rescaling: 27.4%, Approx Optimisation (norm0 = 184.1, norm1 = 225.5, norm01 = 347.8)
#> 2026-07-27 22:32:07.764371: iinla: |lin1-lin0| = 347.8
#>   <eta-lin1,delta>/|delta| = -198.3
#>   |eta-lin0 - delta <delta,eta-lin0>/<delta,delta>| = 107.4
#> 2026-07-27 22:32:07.765726: iinla: Step rescaling: 27.4% (norm0 = 184.1, norm1 = 225.5, norm01 = 347.8)
#> 2026-07-27 22:32:07.767306: iinla: Evaluate predictor linearisation
#> 2026-07-27 22:32:07.769119: Linearise predictor for each observation model
#> 2026-07-27 22:32:07.785129: iinla: Iteration 2 [max: 10]
#> 2026-07-27 22:32:08.534786: iinla: Step rescaling: 162%, Expand (norm0 = 365, norm1 = 141.6, norm01 = 225.4)
#> 2026-07-27 22:32:08.537289: iinla: Step rescaling: 100%, Overstep (norm0 = 225.3, norm1 = 10.71, norm01 = 225.4)
#> 2026-07-27 22:32:08.540334: iinla: Step rescaling: 99.69%, Approx Optimisation (norm0 = 224.6, norm1 = 10.68, norm01 = 225.4)
#> 2026-07-27 22:32:08.54176: iinla: |lin1-lin0| = 225.4
#>   <eta-lin1,delta>/|delta| = -1.007
#>   |eta-lin0 - delta <delta,eta-lin0>/<delta,delta>| = 10.63
#> 2026-07-27 22:32:08.543089: iinla: Step rescaling: 99.7% (norm0 = 224.6, norm1 = 10.68, norm01 = 225.4)
#> 2026-07-27 22:32:08.544678: iinla: Evaluate predictor linearisation
#> 2026-07-27 22:32:08.546407: Linearise predictor for each observation model
#> 2026-07-27 22:32:08.561785: iinla: Max deviation from previous: 46300% of SD, and line search is active
#> [stop if: < 10% and line search inactive]
#> 2026-07-27 22:32:08.563622: iinla: Iteration 3 [max: 10]
#> 2026-07-27 22:32:08.979734: iinla: Step rescaling: 162%, Expand (norm0 = 16.84, norm1 = 6.151, norm01 = 10.68)
#> 2026-07-27 22:32:08.982176: iinla: Step rescaling: 100%, Overstep (norm0 = 10.51, norm1 = 0.1741, norm01 = 10.68)
#> 2026-07-27 22:32:08.985345: iinla: Step rescaling: 101.7%, Approx Optimisation (norm0 = 10.68, norm1 = 0.01176, norm01 = 10.68)
#> 2026-07-27 22:32:08.986944: iinla: |lin1-lin0| = 10.68
#>   <eta-lin1,delta>/|delta| = -1.242e-05
#>   |eta-lin0 - delta <delta,eta-lin0>/<delta,delta>| = 0.01176
#> 2026-07-27 22:32:08.98853: iinla: Step rescaling: 102% (norm0 = 10.68, norm1 = 0.01176, norm01 = 10.68)
#> 2026-07-27 22:32:08.990163: iinla: Evaluate predictor linearisation
#> 2026-07-27 22:32:08.991983: Linearise predictor for each observation model
#> 2026-07-27 22:32:09.007997: iinla: Max deviation from previous: 496% of SD, and line search is active
#> [stop if: < 10% and line search inactive]
#> 2026-07-27 22:32:09.009904: iinla: Iteration 4 [max: 10]
#> 2026-07-27 22:32:09.447935: iinla: Step rescaling: 162%, Expand (norm0 = 0.01902, norm1 = 0.007265, norm01 = 0.01176)
#> 2026-07-27 22:32:09.450835: iinla: Step rescaling: 100%, Overstep (norm0 = 0.01176, norm1 = 3.793e-08, norm01 = 0.01176)
#> 2026-07-27 22:32:09.454018: iinla: Step rescaling: 100%, Approx Optimisation (norm0 = 0.01176, norm1 = 3.784e-08, norm01 = 0.01176)
#> 2026-07-27 22:32:09.455495: iinla: |lin1-lin0| = 0.01176
#>   <eta-lin1,delta>/|delta| = 5.35e-11
#>   |eta-lin0 - delta <delta,eta-lin0>/<delta,delta>| = 3.784e-08
#> 2026-07-27 22:32:09.457227: iinla: Evaluate predictor linearisation
#> 2026-07-27 22:32:09.459122: Linearise predictor for each observation model
#> 2026-07-27 22:32:09.475161: iinla: Max deviation from previous: 8.05% of SD, and line search is inactive
#> [stop if: < 10% and line search inactive]
#> 2026-07-27 22:32:09.476402: iinla: Convergence criterion met.
#>        Running final INLA integration step with known theta mode.
#> 2026-07-27 22:32:09.478259: iinla: Iteration 5 [max: 10]
#> 2026-07-27 22:32:09.913486: iinla: Computation completed
#> 2026-07-27 22:32:11.836655: bru: Preprocessing
#> 2026-07-27 22:32:11.84297: Evaluate component inputs for each observation model
#> 2026-07-27 22:32:11.844427: bru_input(bru_comp_list)
#> 2026-07-27 22:32:11.845762: bru_input.bru_comp(Intercept)
#> 2026-07-27 22:32:11.847157: bru_input.bm_pipe(Intercept)
#> 2026-07-27 22:32:11.848596: bru_input.bm_multi(Intercept:core)
#> 2026-07-27 22:32:11.850006: bru_input.bru_mapper(Intercept:core:main)
#> 2026-07-27 22:32:11.851302: bru_input(bru_input) for (Intercept)
#> 2026-07-27 22:32:11.864159: iinla: Start
#> 2026-07-27 22:32:11.865752: iinla: Evaluate component simplifications
#> 2026-07-27 22:32:11.866909: Simplify component mappers for each observation model
#> 2026-07-27 22:32:11.868214: Simplify component 'NA'
#> 2026-07-27 22:32:11.871889: iinla: Evaluate predictor linearisation
#> 2026-07-27 22:32:11.873771: Linearise predictor for each observation model
#> 2026-07-27 22:32:11.883599: iinla: Construct inla stack
#> 2026-07-27 22:32:11.894137: iinla: Model initialisation completed
#> 2026-07-27 22:32:11.895514: iinla: Iteration 1 [max: 1]
#> 2026-07-27 22:32:12.328147: iinla: Computation completed
#> 2026-07-27 22:32:12.455093: bru: Preprocessing
#> 2026-07-27 22:32:12.463734: Evaluate component inputs for each observation model
#> 2026-07-27 22:32:12.465137: bru_input(bru_comp_list)
#> 2026-07-27 22:32:12.466428: bru_input.bru_comp(Intercept)
#> 2026-07-27 22:32:12.46769: bru_input.bm_pipe(Intercept)
#> 2026-07-27 22:32:12.469057: bru_input.bm_multi(Intercept:core)
#> 2026-07-27 22:32:12.470451: bru_input.bru_mapper(Intercept:core:main)
#> 2026-07-27 22:32:12.471742: bru_input(bru_input) for (Intercept)
#> 2026-07-27 22:32:12.478606: bru_input.bru_comp(field)
#> 2026-07-27 22:32:12.479922: bru_input.bm_pipe(field)
#> 2026-07-27 22:32:12.481318: bru_input.bm_multi(field:core)
#> 2026-07-27 22:32:12.482729: bru_input.bru_mapper(field:core:main)
#> 2026-07-27 22:32:12.483964: bru_input(bru_input) for (field)
#> 2026-07-27 22:32:12.497068: iinla: Start
#> 2026-07-27 22:32:12.498599: iinla: Evaluate component simplifications
#> 2026-07-27 22:32:12.499698: Simplify component mappers for each observation model
#> 2026-07-27 22:32:12.501082: Simplify component 'NA'
#> 2026-07-27 22:32:12.504176: Simplify component 'NA'
#> 2026-07-27 22:32:12.514003: iinla: Evaluate predictor linearisation
#> 2026-07-27 22:32:12.515951: Linearise predictor for each observation model
#> 2026-07-27 22:32:12.523164: iinla: Construct inla stack
#> 2026-07-27 22:32:12.53781: iinla: Model initialisation completed
#> 2026-07-27 22:32:12.53919: iinla: Iteration 1 [max: 1]
#> 2026-07-27 22:32:13.80245: iinla: Computation completed
#> 2026-07-27 22:32:14.953495: bru: Preprocessing
#> 2026-07-27 22:32:14.964985: Evaluate component inputs for each observation model
#> 2026-07-27 22:32:14.96635: bru_input(bru_comp_list)
#> 2026-07-27 22:32:14.96764: bru_input.bru_comp(x)
#> 2026-07-27 22:32:14.968875: bru_input.bm_pipe(x)
#> 2026-07-27 22:32:14.970247: bru_input.bm_multi(x:core)
#> 2026-07-27 22:32:14.971619: bru_input.bru_mapper(x:core:main)
#> 2026-07-27 22:32:14.972856: bru_input(bru_input) for (x)
#> 2026-07-27 22:32:14.97977: bru_input.bru_comp(field)
#> 2026-07-27 22:32:14.981063: bru_input.bm_pipe(field)
#> 2026-07-27 22:32:14.982449: bru_input.bm_multi(field:core)
#> 2026-07-27 22:32:14.983814: bru_input.bru_mapper(field:core:main)
#> 2026-07-27 22:32:14.985028: bru_input(bru_input) for (field)
#> 2026-07-27 22:32:14.998902: iinla: Start
#> 2026-07-27 22:32:15.000453: iinla: Evaluate component simplifications
#> 2026-07-27 22:32:15.001601: Simplify component mappers for each observation model
#> 2026-07-27 22:32:15.002911: Simplify component 'NA'
#> 2026-07-27 22:32:15.005994: Simplify component 'NA'
#> 2026-07-27 22:32:15.009884: iinla: Evaluate predictor linearisation
#> 2026-07-27 22:32:15.011745: Linearise predictor for each observation model
#> 2026-07-27 22:32:15.018406: iinla: Construct inla stack
#> 2026-07-27 22:32:15.031246: iinla: Model initialisation completed
#> 2026-07-27 22:32:15.032713: iinla: Iteration 1 [max: 1]
#> 2026-07-27 22:32:15.464952: iinla: Computation completed
#> 2026-07-27 22:32:17.570656: bru_input(bru_comp_list)
#> 2026-07-27 22:32:19.009765: Evaluate component inputs for each observation model
#> 2026-07-27 22:32:19.010847: bru_input(bru_comp_list)
#> 2026-07-27 22:32:19.01183: bru_input.bru_comp(x)
#> 2026-07-27 22:32:19.012803: bru_input.bm_pipe(x)
#> 2026-07-27 22:32:19.013897: bru_input.bm_multi(x:core)
#> 2026-07-27 22:32:19.014938: bru_input.bru_mapper(x:core:main)
#> 2026-07-27 22:32:19.015873: bru_input(bru_input) for (x)
#> 2026-07-27 22:32:19.022463: bru_input.bru_comp(field)
#> 2026-07-27 22:32:19.023462: bru_input.bm_pipe(field)
#> 2026-07-27 22:32:19.024553: bru_input.bm_multi(field:core)
#> 2026-07-27 22:32:19.0256: bru_input.bru_mapper(field:core:main)
#> 2026-07-27 22:32:19.026495: bru_input(bru_input) for (field)
#> 2026-07-27 22:32:19.031764: bru_input.bru_comp(Intercept)
#> 2026-07-27 22:32:19.032772: bru_input.bm_pipe(Intercept)
#> 2026-07-27 22:32:19.033833: bru_input.bm_multi(Intercept:core)
#> 2026-07-27 22:32:19.034889: bru_input.bru_mapper(Intercept:core:main)
#> 2026-07-27 22:32:19.035818: bru_input(bru_input) for (Intercept)
#> 2026-07-27 22:32:19.0472: bru_input(bru_comp_list)
#> 2026-07-27 22:32:19.04819: bru_input.bru_comp(x)
#> 2026-07-27 22:32:19.049124: bru_input.bm_pipe(x)
#> 2026-07-27 22:32:19.050171: bru_input.bm_multi(x:core)
#> 2026-07-27 22:32:19.05119: bru_input.bru_mapper(x:core:main)
#> 2026-07-27 22:32:19.052086: bru_input(bru_input) for (x)
#> 2026-07-27 22:32:19.057762: bru_input.bru_comp(field)
#> 2026-07-27 22:32:19.058745: bru_input.bm_pipe(field)
#> 2026-07-27 22:32:19.059793: bru_input.bm_multi(field:core)
#> 2026-07-27 22:32:19.060826: bru_input.bru_mapper(field:core:main)
#> 2026-07-27 22:32:19.061747: bru_input(bru_input) for (field)
#> 2026-07-27 22:32:19.067079: bru_input.bru_comp(Intercept)
#> 2026-07-27 22:32:19.06805: bru_input.bm_pipe(Intercept)
#> 2026-07-27 22:32:19.069089: bru_input.bm_multi(Intercept:core)
#> 2026-07-27 22:32:19.07013: bru_input.bru_mapper(Intercept:core:main)
#> 2026-07-27 22:32:19.071042: bru_input(bru_input) for (Intercept)
#> 2026-07-27 22:32:19.076368: Simplify component 'NA'
#> 2026-07-27 22:32:19.079232: Simplify component 'NA'
#> 2026-07-27 22:32:19.08224: Simplify component 'NA'
#> 2026-07-27 22:32:19.085423: bru_input(bru_comp_list)
#> 2026-07-27 22:32:19.086432: bru_input.bru_comp(x)
#> 2026-07-27 22:32:19.08738: bru_input.bm_pipe(x)
#> 2026-07-27 22:32:19.088436: bru_input.bm_multi(x:core)
#> 2026-07-27 22:32:19.08948: bru_input.bru_mapper(x:core:main)
#> 2026-07-27 22:32:19.090393: bru_input(bru_input) for (x)
#> 2026-07-27 22:32:19.096122: bru_input.bru_comp(field)
#> 2026-07-27 22:32:19.097161: bru_input.bm_pipe(field)
#> 2026-07-27 22:32:19.098226: bru_input.bm_multi(field:core)
#> 2026-07-27 22:32:19.099248: bru_input.bru_mapper(field:core:main)
#> 2026-07-27 22:32:19.100171: bru_input(bru_input) for (field)
#> 2026-07-27 22:32:19.105451: bru_input.bru_comp(Intercept)
#> 2026-07-27 22:32:19.106439: bru_input.bm_pipe(Intercept)
#> 2026-07-27 22:32:19.107528: bru_input.bm_multi(Intercept:core)
#> 2026-07-27 22:32:19.108605: bru_input.bru_mapper(Intercept:core:main)
#> 2026-07-27 22:32:19.109528: bru_input(bru_input) for (Intercept)
#> 2026-07-27 22:32:19.114878: Simplify component 'NA'
#> 2026-07-27 22:32:19.117663: Simplify component 'NA'
#> 2026-07-27 22:32:19.120672: Simplify component 'NA'
#> 2026-07-27 22:32:19.125011: bru_input(bru_comp_list)
#> 2026-07-27 22:32:19.12603: bru_input.bru_comp(x)
#> 2026-07-27 22:32:19.126996: bru_input.bm_pipe(x)
#> 2026-07-27 22:32:19.128036: bru_input.bm_multi(x:core)
#> 2026-07-27 22:32:19.129066: bru_input.bru_mapper(x:core:main)
#> 2026-07-27 22:32:19.129991: bru_input(bru_input) for (x)
#> 2026-07-27 22:32:19.135701: bru_input.bru_comp(field)
#> 2026-07-27 22:32:19.136734: bru_input.bm_pipe(field)
#> 2026-07-27 22:32:19.137822: bru_input.bm_multi(field:core)
#> 2026-07-27 22:32:19.138868: bru_input.bru_mapper(field:core:main)
#> 2026-07-27 22:32:19.139808: bru_input(bru_input) for (field)
#> 2026-07-27 22:32:19.145206: bru_input.bru_comp(Intercept)
#> 2026-07-27 22:32:19.146236: bru_input.bm_pipe(Intercept)
#> 2026-07-27 22:32:19.147293: bru_input.bm_multi(Intercept:core)
#> 2026-07-27 22:32:19.148319: bru_input.bru_mapper(Intercept:core:main)
#> 2026-07-27 22:32:19.149237: bru_input(bru_input) for (Intercept)
#> 2026-07-27 22:32:19.154591: Simplify component 'NA'
#> 2026-07-27 22:32:19.157407: Simplify component 'NA'
#> 2026-07-27 22:32:19.160405: Simplify component 'NA'
#> 2026-07-27 22:32:23.209219: bru: Preprocessing
#> 2026-07-27 22:32:23.215774: Evaluate component inputs for each observation model
#> 2026-07-27 22:32:23.217264: bru_input(bru_comp_list)
#> 2026-07-27 22:32:23.218594: bru_input.bru_comp(Intercept)
#> 2026-07-27 22:32:23.219868: bru_input.bm_pipe(Intercept)
#> 2026-07-27 22:32:23.221266: bru_input.bm_multi(Intercept:core)
#> 2026-07-27 22:32:23.222789: bru_input.bru_mapper(Intercept:core:main)
#> 2026-07-27 22:32:23.224045: bru_input(bru_input) for (Intercept)
#> 2026-07-27 22:32:23.240597: iinla: Start
#> 2026-07-27 22:32:23.242131: iinla: Evaluate component simplifications
#> 2026-07-27 22:32:23.243196: Simplify component mappers for each observation model
#> 2026-07-27 22:32:23.244423: Simplify component 'NA'
#> 2026-07-27 22:32:23.248101: iinla: Evaluate predictor linearisation
#> 2026-07-27 22:32:23.249917: Linearise predictor for each observation model
#> 2026-07-27 22:32:23.255009: iinla: Construct inla stack
#> 2026-07-27 22:32:23.263951: iinla: Model initialisation completed
#> 2026-07-27 22:32:23.265291: iinla: Iteration 1 [max: 1]
#> 2026-07-27 22:32:24.44852: iinla: Computation completed
#> 2026-07-27 22:32:25.921234: bru: Preprocessing
#> 2026-07-27 22:32:25.930303: Evaluate component inputs for each observation model
#> 2026-07-27 22:32:25.931807: bru_input(bru_comp_list)
#> 2026-07-27 22:32:25.933111: bru_input.bru_comp(Intercept)
#> 2026-07-27 22:32:25.934406: bru_input.bm_pipe(Intercept)
#> 2026-07-27 22:32:25.935826: bru_input.bm_multi(Intercept:core)
#> 2026-07-27 22:32:25.937171: bru_input.bru_mapper(Intercept:core:main)
#> 2026-07-27 22:32:25.938377: bru_input(bru_input) for (Intercept)
#> 2026-07-27 22:32:25.945391: bru_input.bru_comp(field)
#> 2026-07-27 22:32:25.946775: bru_input.bm_pipe(field)
#> 2026-07-27 22:32:25.948157: bru_input.bm_multi(field:core)
#> 2026-07-27 22:32:25.949495: bru_input.bru_mapper(field:core:main)
#> 2026-07-27 22:32:25.950752: bru_input(bru_input) for (field)
#> 2026-07-27 22:32:25.967126: iinla: Start
#> 2026-07-27 22:32:25.9689: iinla: Evaluate component simplifications
#> 2026-07-27 22:32:25.970016: Simplify component mappers for each observation model
#> 2026-07-27 22:32:25.971287: Simplify component 'NA'
#> 2026-07-27 22:32:25.974457: Simplify component 'NA'
#> 2026-07-27 22:32:26.101429: iinla: Evaluate predictor linearisation
#> 2026-07-27 22:32:26.103301: Linearise predictor for each observation model
#> 2026-07-27 22:32:26.110609: iinla: Construct inla stack
#> 2026-07-27 22:32:26.128464: iinla: Model initialisation completed
#> 2026-07-27 22:32:26.12987: iinla: Iteration 1 [max: 1]
#> 2026-07-27 22:33:09.303148: iinla: Computation completed
#> 2026-07-27 22:33:09.761333: bru: Preprocessing
#> 2026-07-27 22:33:09.781136: Evaluate component inputs for each observation model
#> 2026-07-27 22:33:09.783093: bru_input(bru_comp_list)
#> 2026-07-27 22:33:09.78441: bru_input.bru_comp(x)
#> 2026-07-27 22:33:09.785728: bru_input.bm_pipe(x)
#> 2026-07-27 22:33:09.787112: bru_input.bm_multi(x:core)
#> 2026-07-27 22:33:09.788436: bru_input.bru_mapper(x:core:main)
#> 2026-07-27 22:33:09.789652: bru_input(bru_input) for (x)
#> 2026-07-27 22:33:09.796217: bru_input(bru_comp_list)
#> 2026-07-27 22:33:09.797516: bru_input.bru_comp(x)
#> 2026-07-27 22:33:09.798747: bru_input.bm_pipe(x)
#> 2026-07-27 22:33:09.800063: bru_input.bm_multi(x:core)
#> 2026-07-27 22:33:09.801365: bru_input.bru_mapper(x:core:main)
#> 2026-07-27 22:33:09.802561: bru_input(bru_input) for (x)
#> 2026-07-27 22:33:13.255054: bru_input(bru_input) for (LABEL)
#> 2026-07-27 22:33:14.803875: bru_input_text(bru_input) for (LABEL)
#> 2026-07-27 22:33:14.815768: bru_input.bru_comp(x)
#> 2026-07-27 22:33:14.816918: bru_input_text.bm_pipe(x)
#> 2026-07-27 22:33:14.818089: bru_input_text.bm_multi(x:core)
#> 2026-07-27 22:33:14.819231: bru_input_text.bru_mapper(x:core:main)
#> 2026-07-27 22:33:14.820221: bru_input_text(bru_input) for (x)

bru_log(verbosity = 2L)
#> 2026-07-27 22:31:52.447286: inlabru loaded
#> 2026-07-27 22:31:52.44774: Clear override options
#> 2026-07-27 22:32:05.535344: bru: Preprocessing
#> 2026-07-27 22:32:05.65072: iinla: Iteration 1 [max: 1]
#> 2026-07-27 22:32:06.623265: bru: Preprocessing
#> 2026-07-27 22:32:06.709365: iinla: Iteration 1 [max: 1]
#> 2026-07-27 22:32:07.205802: bru: Preprocessing
#> 2026-07-27 22:32:07.303859: iinla: Iteration 1 [max: 10]
#> 2026-07-27 22:32:07.765726: iinla: Step rescaling: 27.4% (norm0 = 184.1, norm1 = 225.5, norm01 = 347.8)
#> 2026-07-27 22:32:07.785129: iinla: Iteration 2 [max: 10]
#> 2026-07-27 22:32:08.543089: iinla: Step rescaling: 99.7% (norm0 = 224.6, norm1 = 10.68, norm01 = 225.4)
#> 2026-07-27 22:32:08.561785: iinla: Max deviation from previous: 46300% of SD, and line search is active
#> [stop if: < 10% and line search inactive]
#> 2026-07-27 22:32:08.563622: iinla: Iteration 3 [max: 10]
#> 2026-07-27 22:32:08.98853: iinla: Step rescaling: 102% (norm0 = 10.68, norm1 = 0.01176, norm01 = 10.68)
#> 2026-07-27 22:32:09.007997: iinla: Max deviation from previous: 496% of SD, and line search is active
#> [stop if: < 10% and line search inactive]
#> 2026-07-27 22:32:09.009904: iinla: Iteration 4 [max: 10]
#> 2026-07-27 22:32:09.475161: iinla: Max deviation from previous: 8.05% of SD, and line search is inactive
#> [stop if: < 10% and line search inactive]
#> 2026-07-27 22:32:09.476402: iinla: Convergence criterion met.
#>        Running final INLA integration step with known theta mode.
#> 2026-07-27 22:32:09.478259: iinla: Iteration 5 [max: 10]
#> 2026-07-27 22:32:11.836655: bru: Preprocessing
#> 2026-07-27 22:32:11.895514: iinla: Iteration 1 [max: 1]
#> 2026-07-27 22:32:12.455093: bru: Preprocessing
#> 2026-07-27 22:32:12.53919: iinla: Iteration 1 [max: 1]
#> 2026-07-27 22:32:14.953495: bru: Preprocessing
#> 2026-07-27 22:32:15.032713: iinla: Iteration 1 [max: 1]
#> 2026-07-27 22:32:23.209219: bru: Preprocessing
#> 2026-07-27 22:32:23.265291: iinla: Iteration 1 [max: 1]
#> 2026-07-27 22:32:25.921234: bru: Preprocessing
#> 2026-07-27 22:32:26.12987: iinla: Iteration 1 [max: 1]
#> 2026-07-27 22:33:09.761333: bru: Preprocessing
print(bru_log(), timestamp = TRUE, verbosity = TRUE)
#> 2026-07-27 22:31:52.447286: inlabru loaded (level 1)
#> 2026-07-27 22:31:52.44774: Clear override options (level 1)
#> 2026-07-27 22:32:05.535344: bru: Preprocessing (level 1)
#> 2026-07-27 22:32:05.557014: Evaluate component inputs for each observation model (level 3)
#> 2026-07-27 22:32:05.559441: bru_input(bru_comp_list) (level 4)
#> 2026-07-27 22:32:05.561167: bru_input.bru_comp(x) (level 4)
#> 2026-07-27 22:32:05.562665: bru_input.bm_pipe(x) (level 5)
#> 2026-07-27 22:32:05.564245: bru_input.bm_multi(x:core) (level 5)
#> 2026-07-27 22:32:05.565727: bru_input.bru_mapper(x:core:main) (level 5)
#> 2026-07-27 22:32:05.567357: bru_input(bru_input) for (x) (level 5)
#> 2026-07-27 22:32:05.574951: bru_input.bru_comp(Intercept) (level 4)
#> 2026-07-27 22:32:05.576444: bru_input.bm_pipe(Intercept) (level 5)
#> 2026-07-27 22:32:05.578092: bru_input.bm_multi(Intercept:core) (level 5)
#> 2026-07-27 22:32:05.579476: bru_input.bru_mapper(Intercept:core:main) (level 5)
#> 2026-07-27 22:32:05.580757: bru_input(bru_input) for (Intercept) (level 5)
#> 2026-07-27 22:32:05.596993: iinla: Start (level 3)
#> 2026-07-27 22:32:05.59862: iinla: Evaluate component simplifications (level 3)
#> 2026-07-27 22:32:05.599805: Simplify component mappers for each observation model (level 3)
#> 2026-07-27 22:32:05.601259: Simplify component 'NA' (level 4)
#> 2026-07-27 22:32:05.604469: Simplify component 'NA' (level 4)
#> 2026-07-27 22:32:05.607868: iinla: Evaluate predictor linearisation (level 3)
#> 2026-07-27 22:32:05.60991: Linearise predictor for each observation model (level 3)
#> 2026-07-27 22:32:05.630034: iinla: Construct inla stack (level 3)
#> 2026-07-27 22:32:05.649276: iinla: Model initialisation completed (level 3)
#> 2026-07-27 22:32:05.65072: iinla: Iteration 1 [max: 1] (level 1)
#> 2026-07-27 22:32:06.600036: iinla: Computation completed (level 3)
#> 2026-07-27 22:32:06.623265: bru: Preprocessing (level 1)
#> 2026-07-27 22:32:06.632179: Evaluate component inputs for each observation model (level 3)
#> 2026-07-27 22:32:06.633593: bru_input(bru_comp_list) (level 4)
#> 2026-07-27 22:32:06.634923: bru_input.bru_comp(x) (level 4)
#> 2026-07-27 22:32:06.636266: bru_input.bm_pipe(x) (level 5)
#> 2026-07-27 22:32:06.63777: bru_input.bm_multi(x:core) (level 5)
#> 2026-07-27 22:32:06.639199: bru_input.bru_mapper(x:core:main) (level 5)
#> 2026-07-27 22:32:06.640574: bru_input(bru_input) for (x) (level 5)
#> 2026-07-27 22:32:06.648298: bru_input.bru_comp(Intercept) (level 4)
#> 2026-07-27 22:32:06.649708: bru_input.bm_pipe(Intercept) (level 5)
#> 2026-07-27 22:32:06.651115: bru_input.bm_multi(Intercept:core) (level 5)
#> 2026-07-27 22:32:06.652484: bru_input.bru_mapper(Intercept:core:main) (level 5)
#> 2026-07-27 22:32:06.653812: bru_input(bru_input) for (Intercept) (level 5)
#> 2026-07-27 22:32:06.666989: iinla: Start (level 3)
#> 2026-07-27 22:32:06.668622: iinla: Evaluate component simplifications (level 3)
#> 2026-07-27 22:32:06.669775: Simplify component mappers for each observation model (level 3)
#> 2026-07-27 22:32:06.671082: Simplify component 'NA' (level 4)
#> 2026-07-27 22:32:06.674181: Simplify component 'NA' (level 4)
#> 2026-07-27 22:32:06.678353: iinla: Evaluate predictor linearisation (level 3)
#> 2026-07-27 22:32:06.681037: Linearise predictor for each observation model (level 3)
#> 2026-07-27 22:32:06.690596: iinla: Construct inla stack (level 3)
#> 2026-07-27 22:32:06.707399: iinla: Model initialisation completed (level 3)
#> 2026-07-27 22:32:06.709365: iinla: Iteration 1 [max: 1] (level 1)
#> 2026-07-27 22:32:07.178737: iinla: Computation completed (level 3)
#> 2026-07-27 22:32:07.205802: bru: Preprocessing (level 1)
#> 2026-07-27 22:32:07.216748: Evaluate component inputs for each observation model (level 3)
#> 2026-07-27 22:32:07.218695: bru_input(bru_comp_list) (level 4)
#> 2026-07-27 22:32:07.220516: bru_input.bru_comp(z) (level 4)
#> 2026-07-27 22:32:07.222291: bru_input.bm_pipe(z) (level 5)
#> 2026-07-27 22:32:07.224382: bru_input.bm_multi(z:core) (level 5)
#> 2026-07-27 22:32:07.226321: bru_input.bru_mapper(z:core:main) (level 5)
#> 2026-07-27 22:32:07.228101: bru_input(bru_input) for (z) (level 5)
#> 2026-07-27 22:32:07.235719: bru_input.bru_comp(Intercept) (level 4)
#> 2026-07-27 22:32:07.237532: bru_input.bm_pipe(Intercept) (level 5)
#> 2026-07-27 22:32:07.239467: bru_input.bm_multi(Intercept:core) (level 5)
#> 2026-07-27 22:32:07.241391: bru_input.bru_mapper(Intercept:core:main) (level 5)
#> 2026-07-27 22:32:07.243237: bru_input(bru_input) for (Intercept) (level 5)
#> 2026-07-27 22:32:07.258339: iinla: Start (level 3)
#> 2026-07-27 22:32:07.260482: iinla: Evaluate component simplifications (level 3)
#> 2026-07-27 22:32:07.262089: Simplify component mappers for each observation model (level 3)
#> 2026-07-27 22:32:07.263993: Simplify component 'NA' (level 4)
#> 2026-07-27 22:32:07.26787: Simplify component 'NA' (level 4)
#> 2026-07-27 22:32:07.272138: iinla: Evaluate predictor linearisation (level 3)
#> 2026-07-27 22:32:07.274735: Linearise predictor for each observation model (level 3)
#> 2026-07-27 22:32:07.284394: iinla: Construct inla stack (level 3)
#> 2026-07-27 22:32:07.301619: iinla: Model initialisation completed (level 3)
#> 2026-07-27 22:32:07.303859: iinla: Iteration 1 [max: 10] (level 1)
#> 2026-07-27 22:32:07.757508: iinla: Step rescaling: 61.8%, Contract (norm0 = 1969, norm1 = 1807, norm01 = 347.8) (level 3)
#> 2026-07-27 22:32:07.759907: iinla: Step rescaling: 38.2%, Contract (norm0 = 398.9, norm1 = 295.1, norm01 = 347.8) (level 3)
#> 2026-07-27 22:32:07.762973: iinla: Step rescaling: 27.4%, Approx Optimisation (norm0 = 184.1, norm1 = 225.5, norm01 = 347.8) (level 3)
#> 2026-07-27 22:32:07.764371: iinla: |lin1-lin0| = 347.8
#>   <eta-lin1,delta>/|delta| = -198.3
#>   |eta-lin0 - delta <delta,eta-lin0>/<delta,delta>| = 107.4 (level 4)
#> 2026-07-27 22:32:07.765726: iinla: Step rescaling: 27.4% (norm0 = 184.1, norm1 = 225.5, norm01 = 347.8) (level 2)
#> 2026-07-27 22:32:07.767306: iinla: Evaluate predictor linearisation (level 3)
#> 2026-07-27 22:32:07.769119: Linearise predictor for each observation model (level 3)
#> 2026-07-27 22:32:07.785129: iinla: Iteration 2 [max: 10] (level 1)
#> 2026-07-27 22:32:08.534786: iinla: Step rescaling: 162%, Expand (norm0 = 365, norm1 = 141.6, norm01 = 225.4) (level 3)
#> 2026-07-27 22:32:08.537289: iinla: Step rescaling: 100%, Overstep (norm0 = 225.3, norm1 = 10.71, norm01 = 225.4) (level 3)
#> 2026-07-27 22:32:08.540334: iinla: Step rescaling: 99.69%, Approx Optimisation (norm0 = 224.6, norm1 = 10.68, norm01 = 225.4) (level 3)
#> 2026-07-27 22:32:08.54176: iinla: |lin1-lin0| = 225.4
#>   <eta-lin1,delta>/|delta| = -1.007
#>   |eta-lin0 - delta <delta,eta-lin0>/<delta,delta>| = 10.63 (level 4)
#> 2026-07-27 22:32:08.543089: iinla: Step rescaling: 99.7% (norm0 = 224.6, norm1 = 10.68, norm01 = 225.4) (level 2)
#> 2026-07-27 22:32:08.544678: iinla: Evaluate predictor linearisation (level 3)
#> 2026-07-27 22:32:08.546407: Linearise predictor for each observation model (level 3)
#> 2026-07-27 22:32:08.561785: iinla: Max deviation from previous: 46300% of SD, and line search is active
#> [stop if: < 10% and line search inactive] (level 1)
#> 2026-07-27 22:32:08.563622: iinla: Iteration 3 [max: 10] (level 1)
#> 2026-07-27 22:32:08.979734: iinla: Step rescaling: 162%, Expand (norm0 = 16.84, norm1 = 6.151, norm01 = 10.68) (level 3)
#> 2026-07-27 22:32:08.982176: iinla: Step rescaling: 100%, Overstep (norm0 = 10.51, norm1 = 0.1741, norm01 = 10.68) (level 3)
#> 2026-07-27 22:32:08.985345: iinla: Step rescaling: 101.7%, Approx Optimisation (norm0 = 10.68, norm1 = 0.01176, norm01 = 10.68) (level 3)
#> 2026-07-27 22:32:08.986944: iinla: |lin1-lin0| = 10.68
#>   <eta-lin1,delta>/|delta| = -1.242e-05
#>   |eta-lin0 - delta <delta,eta-lin0>/<delta,delta>| = 0.01176 (level 4)
#> 2026-07-27 22:32:08.98853: iinla: Step rescaling: 102% (norm0 = 10.68, norm1 = 0.01176, norm01 = 10.68) (level 2)
#> 2026-07-27 22:32:08.990163: iinla: Evaluate predictor linearisation (level 3)
#> 2026-07-27 22:32:08.991983: Linearise predictor for each observation model (level 3)
#> 2026-07-27 22:32:09.007997: iinla: Max deviation from previous: 496% of SD, and line search is active
#> [stop if: < 10% and line search inactive] (level 1)
#> 2026-07-27 22:32:09.009904: iinla: Iteration 4 [max: 10] (level 1)
#> 2026-07-27 22:32:09.447935: iinla: Step rescaling: 162%, Expand (norm0 = 0.01902, norm1 = 0.007265, norm01 = 0.01176) (level 3)
#> 2026-07-27 22:32:09.450835: iinla: Step rescaling: 100%, Overstep (norm0 = 0.01176, norm1 = 3.793e-08, norm01 = 0.01176) (level 3)
#> 2026-07-27 22:32:09.454018: iinla: Step rescaling: 100%, Approx Optimisation (norm0 = 0.01176, norm1 = 3.784e-08, norm01 = 0.01176) (level 3)
#> 2026-07-27 22:32:09.455495: iinla: |lin1-lin0| = 0.01176
#>   <eta-lin1,delta>/|delta| = 5.35e-11
#>   |eta-lin0 - delta <delta,eta-lin0>/<delta,delta>| = 3.784e-08 (level 4)
#> 2026-07-27 22:32:09.457227: iinla: Evaluate predictor linearisation (level 3)
#> 2026-07-27 22:32:09.459122: Linearise predictor for each observation model (level 3)
#> 2026-07-27 22:32:09.475161: iinla: Max deviation from previous: 8.05% of SD, and line search is inactive
#> [stop if: < 10% and line search inactive] (level 1)
#> 2026-07-27 22:32:09.476402: iinla: Convergence criterion met.
#>        Running final INLA integration step with known theta mode. (level 1)
#> 2026-07-27 22:32:09.478259: iinla: Iteration 5 [max: 10] (level 1)
#> 2026-07-27 22:32:09.913486: iinla: Computation completed (level 3)
#> 2026-07-27 22:32:11.836655: bru: Preprocessing (level 1)
#> 2026-07-27 22:32:11.84297: Evaluate component inputs for each observation model (level 3)
#> 2026-07-27 22:32:11.844427: bru_input(bru_comp_list) (level 4)
#> 2026-07-27 22:32:11.845762: bru_input.bru_comp(Intercept) (level 4)
#> 2026-07-27 22:32:11.847157: bru_input.bm_pipe(Intercept) (level 5)
#> 2026-07-27 22:32:11.848596: bru_input.bm_multi(Intercept:core) (level 5)
#> 2026-07-27 22:32:11.850006: bru_input.bru_mapper(Intercept:core:main) (level 5)
#> 2026-07-27 22:32:11.851302: bru_input(bru_input) for (Intercept) (level 5)
#> 2026-07-27 22:32:11.864159: iinla: Start (level 3)
#> 2026-07-27 22:32:11.865752: iinla: Evaluate component simplifications (level 3)
#> 2026-07-27 22:32:11.866909: Simplify component mappers for each observation model (level 3)
#> 2026-07-27 22:32:11.868214: Simplify component 'NA' (level 4)
#> 2026-07-27 22:32:11.871889: iinla: Evaluate predictor linearisation (level 3)
#> 2026-07-27 22:32:11.873771: Linearise predictor for each observation model (level 3)
#> 2026-07-27 22:32:11.883599: iinla: Construct inla stack (level 3)
#> 2026-07-27 22:32:11.894137: iinla: Model initialisation completed (level 3)
#> 2026-07-27 22:32:11.895514: iinla: Iteration 1 [max: 1] (level 1)
#> 2026-07-27 22:32:12.328147: iinla: Computation completed (level 3)
#> 2026-07-27 22:32:12.455093: bru: Preprocessing (level 1)
#> 2026-07-27 22:32:12.463734: Evaluate component inputs for each observation model (level 3)
#> 2026-07-27 22:32:12.465137: bru_input(bru_comp_list) (level 4)
#> 2026-07-27 22:32:12.466428: bru_input.bru_comp(Intercept) (level 4)
#> 2026-07-27 22:32:12.46769: bru_input.bm_pipe(Intercept) (level 5)
#> 2026-07-27 22:32:12.469057: bru_input.bm_multi(Intercept:core) (level 5)
#> 2026-07-27 22:32:12.470451: bru_input.bru_mapper(Intercept:core:main) (level 5)
#> 2026-07-27 22:32:12.471742: bru_input(bru_input) for (Intercept) (level 5)
#> 2026-07-27 22:32:12.478606: bru_input.bru_comp(field) (level 4)
#> 2026-07-27 22:32:12.479922: bru_input.bm_pipe(field) (level 5)
#> 2026-07-27 22:32:12.481318: bru_input.bm_multi(field:core) (level 5)
#> 2026-07-27 22:32:12.482729: bru_input.bru_mapper(field:core:main) (level 5)
#> 2026-07-27 22:32:12.483964: bru_input(bru_input) for (field) (level 5)
#> 2026-07-27 22:32:12.497068: iinla: Start (level 3)
#> 2026-07-27 22:32:12.498599: iinla: Evaluate component simplifications (level 3)
#> 2026-07-27 22:32:12.499698: Simplify component mappers for each observation model (level 3)
#> 2026-07-27 22:32:12.501082: Simplify component 'NA' (level 4)
#> 2026-07-27 22:32:12.504176: Simplify component 'NA' (level 4)
#> 2026-07-27 22:32:12.514003: iinla: Evaluate predictor linearisation (level 3)
#> 2026-07-27 22:32:12.515951: Linearise predictor for each observation model (level 3)
#> 2026-07-27 22:32:12.523164: iinla: Construct inla stack (level 3)
#> 2026-07-27 22:32:12.53781: iinla: Model initialisation completed (level 3)
#> 2026-07-27 22:32:12.53919: iinla: Iteration 1 [max: 1] (level 1)
#> 2026-07-27 22:32:13.80245: iinla: Computation completed (level 3)
#> 2026-07-27 22:32:14.953495: bru: Preprocessing (level 1)
#> 2026-07-27 22:32:14.964985: Evaluate component inputs for each observation model (level 3)
#> 2026-07-27 22:32:14.96635: bru_input(bru_comp_list) (level 4)
#> 2026-07-27 22:32:14.96764: bru_input.bru_comp(x) (level 4)
#> 2026-07-27 22:32:14.968875: bru_input.bm_pipe(x) (level 5)
#> 2026-07-27 22:32:14.970247: bru_input.bm_multi(x:core) (level 5)
#> 2026-07-27 22:32:14.971619: bru_input.bru_mapper(x:core:main) (level 5)
#> 2026-07-27 22:32:14.972856: bru_input(bru_input) for (x) (level 5)
#> 2026-07-27 22:32:14.97977: bru_input.bru_comp(field) (level 4)
#> 2026-07-27 22:32:14.981063: bru_input.bm_pipe(field) (level 5)
#> 2026-07-27 22:32:14.982449: bru_input.bm_multi(field:core) (level 5)
#> 2026-07-27 22:32:14.983814: bru_input.bru_mapper(field:core:main) (level 5)
#> 2026-07-27 22:32:14.985028: bru_input(bru_input) for (field) (level 5)
#> 2026-07-27 22:32:14.998902: iinla: Start (level 3)
#> 2026-07-27 22:32:15.000453: iinla: Evaluate component simplifications (level 3)
#> 2026-07-27 22:32:15.001601: Simplify component mappers for each observation model (level 3)
#> 2026-07-27 22:32:15.002911: Simplify component 'NA' (level 4)
#> 2026-07-27 22:32:15.005994: Simplify component 'NA' (level 4)
#> 2026-07-27 22:32:15.009884: iinla: Evaluate predictor linearisation (level 3)
#> 2026-07-27 22:32:15.011745: Linearise predictor for each observation model (level 3)
#> 2026-07-27 22:32:15.018406: iinla: Construct inla stack (level 3)
#> 2026-07-27 22:32:15.031246: iinla: Model initialisation completed (level 3)
#> 2026-07-27 22:32:15.032713: iinla: Iteration 1 [max: 1] (level 1)
#> 2026-07-27 22:32:15.464952: iinla: Computation completed (level 3)
#> 2026-07-27 22:32:17.570656: bru_input(bru_comp_list) (level 4)
#> 2026-07-27 22:32:19.009765: Evaluate component inputs for each observation model (level 3)
#> 2026-07-27 22:32:19.010847: bru_input(bru_comp_list) (level 4)
#> 2026-07-27 22:32:19.01183: bru_input.bru_comp(x) (level 4)
#> 2026-07-27 22:32:19.012803: bru_input.bm_pipe(x) (level 5)
#> 2026-07-27 22:32:19.013897: bru_input.bm_multi(x:core) (level 5)
#> 2026-07-27 22:32:19.014938: bru_input.bru_mapper(x:core:main) (level 5)
#> 2026-07-27 22:32:19.015873: bru_input(bru_input) for (x) (level 5)
#> 2026-07-27 22:32:19.022463: bru_input.bru_comp(field) (level 4)
#> 2026-07-27 22:32:19.023462: bru_input.bm_pipe(field) (level 5)
#> 2026-07-27 22:32:19.024553: bru_input.bm_multi(field:core) (level 5)
#> 2026-07-27 22:32:19.0256: bru_input.bru_mapper(field:core:main) (level 5)
#> 2026-07-27 22:32:19.026495: bru_input(bru_input) for (field) (level 5)
#> 2026-07-27 22:32:19.031764: bru_input.bru_comp(Intercept) (level 4)
#> 2026-07-27 22:32:19.032772: bru_input.bm_pipe(Intercept) (level 5)
#> 2026-07-27 22:32:19.033833: bru_input.bm_multi(Intercept:core) (level 5)
#> 2026-07-27 22:32:19.034889: bru_input.bru_mapper(Intercept:core:main) (level 5)
#> 2026-07-27 22:32:19.035818: bru_input(bru_input) for (Intercept) (level 5)
#> 2026-07-27 22:32:19.0472: bru_input(bru_comp_list) (level 4)
#> 2026-07-27 22:32:19.04819: bru_input.bru_comp(x) (level 4)
#> 2026-07-27 22:32:19.049124: bru_input.bm_pipe(x) (level 5)
#> 2026-07-27 22:32:19.050171: bru_input.bm_multi(x:core) (level 5)
#> 2026-07-27 22:32:19.05119: bru_input.bru_mapper(x:core:main) (level 5)
#> 2026-07-27 22:32:19.052086: bru_input(bru_input) for (x) (level 5)
#> 2026-07-27 22:32:19.057762: bru_input.bru_comp(field) (level 4)
#> 2026-07-27 22:32:19.058745: bru_input.bm_pipe(field) (level 5)
#> 2026-07-27 22:32:19.059793: bru_input.bm_multi(field:core) (level 5)
#> 2026-07-27 22:32:19.060826: bru_input.bru_mapper(field:core:main) (level 5)
#> 2026-07-27 22:32:19.061747: bru_input(bru_input) for (field) (level 5)
#> 2026-07-27 22:32:19.067079: bru_input.bru_comp(Intercept) (level 4)
#> 2026-07-27 22:32:19.06805: bru_input.bm_pipe(Intercept) (level 5)
#> 2026-07-27 22:32:19.069089: bru_input.bm_multi(Intercept:core) (level 5)
#> 2026-07-27 22:32:19.07013: bru_input.bru_mapper(Intercept:core:main) (level 5)
#> 2026-07-27 22:32:19.071042: bru_input(bru_input) for (Intercept) (level 5)
#> 2026-07-27 22:32:19.076368: Simplify component 'NA' (level 4)
#> 2026-07-27 22:32:19.079232: Simplify component 'NA' (level 4)
#> 2026-07-27 22:32:19.08224: Simplify component 'NA' (level 4)
#> 2026-07-27 22:32:19.085423: bru_input(bru_comp_list) (level 4)
#> 2026-07-27 22:32:19.086432: bru_input.bru_comp(x) (level 4)
#> 2026-07-27 22:32:19.08738: bru_input.bm_pipe(x) (level 5)
#> 2026-07-27 22:32:19.088436: bru_input.bm_multi(x:core) (level 5)
#> 2026-07-27 22:32:19.08948: bru_input.bru_mapper(x:core:main) (level 5)
#> 2026-07-27 22:32:19.090393: bru_input(bru_input) for (x) (level 5)
#> 2026-07-27 22:32:19.096122: bru_input.bru_comp(field) (level 4)
#> 2026-07-27 22:32:19.097161: bru_input.bm_pipe(field) (level 5)
#> 2026-07-27 22:32:19.098226: bru_input.bm_multi(field:core) (level 5)
#> 2026-07-27 22:32:19.099248: bru_input.bru_mapper(field:core:main) (level 5)
#> 2026-07-27 22:32:19.100171: bru_input(bru_input) for (field) (level 5)
#> 2026-07-27 22:32:19.105451: bru_input.bru_comp(Intercept) (level 4)
#> 2026-07-27 22:32:19.106439: bru_input.bm_pipe(Intercept) (level 5)
#> 2026-07-27 22:32:19.107528: bru_input.bm_multi(Intercept:core) (level 5)
#> 2026-07-27 22:32:19.108605: bru_input.bru_mapper(Intercept:core:main) (level 5)
#> 2026-07-27 22:32:19.109528: bru_input(bru_input) for (Intercept) (level 5)
#> 2026-07-27 22:32:19.114878: Simplify component 'NA' (level 4)
#> 2026-07-27 22:32:19.117663: Simplify component 'NA' (level 4)
#> 2026-07-27 22:32:19.120672: Simplify component 'NA' (level 4)
#> 2026-07-27 22:32:19.125011: bru_input(bru_comp_list) (level 4)
#> 2026-07-27 22:32:19.12603: bru_input.bru_comp(x) (level 4)
#> 2026-07-27 22:32:19.126996: bru_input.bm_pipe(x) (level 5)
#> 2026-07-27 22:32:19.128036: bru_input.bm_multi(x:core) (level 5)
#> 2026-07-27 22:32:19.129066: bru_input.bru_mapper(x:core:main) (level 5)
#> 2026-07-27 22:32:19.129991: bru_input(bru_input) for (x) (level 5)
#> 2026-07-27 22:32:19.135701: bru_input.bru_comp(field) (level 4)
#> 2026-07-27 22:32:19.136734: bru_input.bm_pipe(field) (level 5)
#> 2026-07-27 22:32:19.137822: bru_input.bm_multi(field:core) (level 5)
#> 2026-07-27 22:32:19.138868: bru_input.bru_mapper(field:core:main) (level 5)
#> 2026-07-27 22:32:19.139808: bru_input(bru_input) for (field) (level 5)
#> 2026-07-27 22:32:19.145206: bru_input.bru_comp(Intercept) (level 4)
#> 2026-07-27 22:32:19.146236: bru_input.bm_pipe(Intercept) (level 5)
#> 2026-07-27 22:32:19.147293: bru_input.bm_multi(Intercept:core) (level 5)
#> 2026-07-27 22:32:19.148319: bru_input.bru_mapper(Intercept:core:main) (level 5)
#> 2026-07-27 22:32:19.149237: bru_input(bru_input) for (Intercept) (level 5)
#> 2026-07-27 22:32:19.154591: Simplify component 'NA' (level 4)
#> 2026-07-27 22:32:19.157407: Simplify component 'NA' (level 4)
#> 2026-07-27 22:32:19.160405: Simplify component 'NA' (level 4)
#> 2026-07-27 22:32:23.209219: bru: Preprocessing (level 1)
#> 2026-07-27 22:32:23.215774: Evaluate component inputs for each observation model (level 3)
#> 2026-07-27 22:32:23.217264: bru_input(bru_comp_list) (level 4)
#> 2026-07-27 22:32:23.218594: bru_input.bru_comp(Intercept) (level 4)
#> 2026-07-27 22:32:23.219868: bru_input.bm_pipe(Intercept) (level 5)
#> 2026-07-27 22:32:23.221266: bru_input.bm_multi(Intercept:core) (level 5)
#> 2026-07-27 22:32:23.222789: bru_input.bru_mapper(Intercept:core:main) (level 5)
#> 2026-07-27 22:32:23.224045: bru_input(bru_input) for (Intercept) (level 5)
#> 2026-07-27 22:32:23.240597: iinla: Start (level 3)
#> 2026-07-27 22:32:23.242131: iinla: Evaluate component simplifications (level 3)
#> 2026-07-27 22:32:23.243196: Simplify component mappers for each observation model (level 3)
#> 2026-07-27 22:32:23.244423: Simplify component 'NA' (level 4)
#> 2026-07-27 22:32:23.248101: iinla: Evaluate predictor linearisation (level 3)
#> 2026-07-27 22:32:23.249917: Linearise predictor for each observation model (level 3)
#> 2026-07-27 22:32:23.255009: iinla: Construct inla stack (level 3)
#> 2026-07-27 22:32:23.263951: iinla: Model initialisation completed (level 3)
#> 2026-07-27 22:32:23.265291: iinla: Iteration 1 [max: 1] (level 1)
#> 2026-07-27 22:32:24.44852: iinla: Computation completed (level 3)
#> 2026-07-27 22:32:25.921234: bru: Preprocessing (level 1)
#> 2026-07-27 22:32:25.930303: Evaluate component inputs for each observation model (level 3)
#> 2026-07-27 22:32:25.931807: bru_input(bru_comp_list) (level 4)
#> 2026-07-27 22:32:25.933111: bru_input.bru_comp(Intercept) (level 4)
#> 2026-07-27 22:32:25.934406: bru_input.bm_pipe(Intercept) (level 5)
#> 2026-07-27 22:32:25.935826: bru_input.bm_multi(Intercept:core) (level 5)
#> 2026-07-27 22:32:25.937171: bru_input.bru_mapper(Intercept:core:main) (level 5)
#> 2026-07-27 22:32:25.938377: bru_input(bru_input) for (Intercept) (level 5)
#> 2026-07-27 22:32:25.945391: bru_input.bru_comp(field) (level 4)
#> 2026-07-27 22:32:25.946775: bru_input.bm_pipe(field) (level 5)
#> 2026-07-27 22:32:25.948157: bru_input.bm_multi(field:core) (level 5)
#> 2026-07-27 22:32:25.949495: bru_input.bru_mapper(field:core:main) (level 5)
#> 2026-07-27 22:32:25.950752: bru_input(bru_input) for (field) (level 5)
#> 2026-07-27 22:32:25.967126: iinla: Start (level 3)
#> 2026-07-27 22:32:25.9689: iinla: Evaluate component simplifications (level 3)
#> 2026-07-27 22:32:25.970016: Simplify component mappers for each observation model (level 3)
#> 2026-07-27 22:32:25.971287: Simplify component 'NA' (level 4)
#> 2026-07-27 22:32:25.974457: Simplify component 'NA' (level 4)
#> 2026-07-27 22:32:26.101429: iinla: Evaluate predictor linearisation (level 3)
#> 2026-07-27 22:32:26.103301: Linearise predictor for each observation model (level 3)
#> 2026-07-27 22:32:26.110609: iinla: Construct inla stack (level 3)
#> 2026-07-27 22:32:26.128464: iinla: Model initialisation completed (level 3)
#> 2026-07-27 22:32:26.12987: iinla: Iteration 1 [max: 1] (level 1)
#> 2026-07-27 22:33:09.303148: iinla: Computation completed (level 3)
#> 2026-07-27 22:33:09.761333: bru: Preprocessing (level 1)
#> 2026-07-27 22:33:09.781136: Evaluate component inputs for each observation model (level 3)
#> 2026-07-27 22:33:09.783093: bru_input(bru_comp_list) (level 4)
#> 2026-07-27 22:33:09.78441: bru_input.bru_comp(x) (level 4)
#> 2026-07-27 22:33:09.785728: bru_input.bm_pipe(x) (level 5)
#> 2026-07-27 22:33:09.787112: bru_input.bm_multi(x:core) (level 5)
#> 2026-07-27 22:33:09.788436: bru_input.bru_mapper(x:core:main) (level 5)
#> 2026-07-27 22:33:09.789652: bru_input(bru_input) for (x) (level 5)
#> 2026-07-27 22:33:09.796217: bru_input(bru_comp_list) (level 4)
#> 2026-07-27 22:33:09.797516: bru_input.bru_comp(x) (level 4)
#> 2026-07-27 22:33:09.798747: bru_input.bm_pipe(x) (level 5)
#> 2026-07-27 22:33:09.800063: bru_input.bm_multi(x:core) (level 5)
#> 2026-07-27 22:33:09.801365: bru_input.bru_mapper(x:core:main) (level 5)
#> 2026-07-27 22:33:09.802561: bru_input(bru_input) for (x) (level 5)
#> 2026-07-27 22:33:13.255054: bru_input(bru_input) for (LABEL) (level 5)
#> 2026-07-27 22:33:14.803875: bru_input_text(bru_input) for (LABEL) (level 5)
#> 2026-07-27 22:33:14.815768: bru_input.bru_comp(x) (level 4)
#> 2026-07-27 22:33:14.816918: bru_input_text.bm_pipe(x) (level 5)
#> 2026-07-27 22:33:14.818089: bru_input_text.bm_multi(x:core) (level 5)
#> 2026-07-27 22:33:14.819231: bru_input_text.bru_mapper(x:core:main) (level 5)
#> 2026-07-27 22:33:14.820221: bru_input_text(bru_input) for (x) (level 5)
```
