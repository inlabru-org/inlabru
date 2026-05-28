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
[`bru_log_new()`](https://inlabru-org.github.io/inlabru/reference/bru_log_new.md),
[`bru_log_offset()`](https://inlabru-org.github.io/inlabru/reference/bru_log_offset.md),
[`bru_log_reset()`](https://inlabru-org.github.io/inlabru/reference/bru_log_reset.md)

## Examples

``` r
bru_log(verbosity = 2L)
#> 2026-05-26 12:45:38.190956: inlabru loaded
#> 2026-05-26 12:45:38.191377: Clear override options
#> 2026-05-26 12:45:49.727601: bru: Preprocessing
#> 2026-05-26 12:45:49.84122: iinla: Iteration 1 [max: 1]
#> 2026-05-26 12:45:50.418018: bru: Preprocessing
#> 2026-05-26 12:45:50.504482: iinla: Iteration 1 [max: 1]
#> 2026-05-26 12:45:50.888623: bru: Preprocessing
#> 2026-05-26 12:45:50.981709: iinla: Iteration 1 [max: 10]
#> 2026-05-26 12:45:51.340373: iinla: Step rescaling: 27.4% (norm0 = 184.1, norm1 = 225.5, norm01 = 347.8)
#> 2026-05-26 12:45:51.368339: iinla: Iteration 2 [max: 10]
#> 2026-05-26 12:45:51.95141: iinla: Step rescaling: 99.7% (norm0 = 224.7, norm1 = 10.68, norm01 = 225.4)
#> 2026-05-26 12:45:51.978858: iinla: Max deviation from previous: 46300% of SD, and line search is active
#> [stop if: < 10% and line search inactive]
#> 2026-05-26 12:45:51.981334: iinla: Iteration 3 [max: 10]
#> 2026-05-26 12:45:52.334209: iinla: Step rescaling: 102% (norm0 = 10.69, norm1 = 0.01176, norm01 = 10.69)
#> 2026-05-26 12:45:52.363075: iinla: Max deviation from previous: 496% of SD, and line search is active
#> [stop if: < 10% and line search inactive]
#> 2026-05-26 12:45:52.365527: iinla: Iteration 4 [max: 10]
#> 2026-05-26 12:45:52.773334: iinla: Max deviation from previous: 8.06% of SD, and line search is inactive
#> [stop if: < 10% and line search inactive]
#> 2026-05-26 12:45:52.77474: iinla: Convergence criterion met.
#>        Running final INLA integration step with known theta mode.
#> 2026-05-26 12:45:52.777451: iinla: Iteration 5 [max: 10]
#> 2026-05-26 12:45:54.292735: bru: Preprocessing
#> 2026-05-26 12:45:54.358945: iinla: Iteration 1 [max: 1]
#> 2026-05-26 12:45:59.674772: bru: Preprocessing
format(bru_log())
#> 2026-05-26 12:45:38.190956: inlabru loaded
#> 2026-05-26 12:45:38.191377: Clear override options
#> 2026-05-26 12:45:49.727601: bru: Preprocessing
#> 2026-05-26 12:45:49.743866: Evaluate component inputs for each observation model
#> 2026-05-26 12:45:49.745247: bru_input(bru_comp_list)
#> 2026-05-26 12:45:49.746636: bru_input.bru_comp(x)
#> 2026-05-26 12:45:49.748101: bru_input.bm_pipe(x)
#> 2026-05-26 12:45:49.749836: bru_input.bm_multi(x:core)
#> 2026-05-26 12:45:49.751425: bru_input.bru_mapper(x:core:main)
#> 2026-05-26 12:45:49.753049: bru_input(bru_input) for (x)
#> 2026-05-26 12:45:49.761505: bru_input.bru_comp(Intercept)
#> 2026-05-26 12:45:49.762754: bru_input.bm_pipe(Intercept)
#> 2026-05-26 12:45:49.764135: bru_input.bm_multi(Intercept:core)
#> 2026-05-26 12:45:49.765461: bru_input.bru_mapper(Intercept:core:main)
#> 2026-05-26 12:45:49.766691: bru_input(bru_input) for (Intercept)
#> 2026-05-26 12:45:49.78371: iinla: Start
#> 2026-05-26 12:45:49.785265: iinla: Evaluate component linearisations
#> 2026-05-26 12:45:49.786408: Linearise components for each observation model
#> 2026-05-26 12:45:49.787807: Linearise component 'x'
#> 2026-05-26 12:45:49.790721: Linearise component 'Intercept'
#> 2026-05-26 12:45:49.793379: iinla: Evaluate component simplifications
#> 2026-05-26 12:45:49.794529: Simplify component mappers for each observation model
#> 2026-05-26 12:45:49.795882: Simplify component 'x'
#> 2026-05-26 12:45:49.798847: Simplify component 'Intercept'
#> 2026-05-26 12:45:49.801742: iinla: Evaluate predictor linearisation
#> 2026-05-26 12:45:49.815638: iinla: Construct inla stack
#> 2026-05-26 12:45:49.839587: iinla: Model initialisation completed
#> 2026-05-26 12:45:49.84122: iinla: Iteration 1 [max: 1]
#> 2026-05-26 12:45:50.39276: iinla: Computation completed
#> 2026-05-26 12:45:50.418018: bru: Preprocessing
#> 2026-05-26 12:45:50.426873: Evaluate component inputs for each observation model
#> 2026-05-26 12:45:50.428097: bru_input(bru_comp_list)
#> 2026-05-26 12:45:50.429435: bru_input.bru_comp(x)
#> 2026-05-26 12:45:50.430851: bru_input.bm_pipe(x)
#> 2026-05-26 12:45:50.432316: bru_input.bm_multi(x:core)
#> 2026-05-26 12:45:50.433751: bru_input.bru_mapper(x:core:main)
#> 2026-05-26 12:45:50.435068: bru_input(bru_input) for (x)
#> 2026-05-26 12:45:50.445308: bru_input.bru_comp(Intercept)
#> 2026-05-26 12:45:50.446794: bru_input.bm_pipe(Intercept)
#> 2026-05-26 12:45:50.448344: bru_input.bm_multi(Intercept:core)
#> 2026-05-26 12:45:50.44979: bru_input.bru_mapper(Intercept:core:main)
#> 2026-05-26 12:45:50.451172: bru_input(bru_input) for (Intercept)
#> 2026-05-26 12:45:50.464854: iinla: Start
#> 2026-05-26 12:45:50.466457: iinla: Evaluate component linearisations
#> 2026-05-26 12:45:50.467615: Linearise components for each observation model
#> 2026-05-26 12:45:50.46896: Linearise component 'x'
#> 2026-05-26 12:45:50.472423: Linearise component 'Intercept'
#> 2026-05-26 12:45:50.475768: iinla: Evaluate component simplifications
#> 2026-05-26 12:45:50.477107: Simplify component mappers for each observation model
#> 2026-05-26 12:45:50.47861: Simplify component 'x'
#> 2026-05-26 12:45:50.481707: Simplify component 'Intercept'
#> 2026-05-26 12:45:50.484711: iinla: Evaluate predictor linearisation
#> 2026-05-26 12:45:50.491022: iinla: Construct inla stack
#> 2026-05-26 12:45:50.502785: iinla: Model initialisation completed
#> 2026-05-26 12:45:50.504482: iinla: Iteration 1 [max: 1]
#> 2026-05-26 12:45:50.867193: iinla: Computation completed
#> 2026-05-26 12:45:50.888623: bru: Preprocessing
#> 2026-05-26 12:45:50.896928: Evaluate component inputs for each observation model
#> 2026-05-26 12:45:50.898171: bru_input(bru_comp_list)
#> 2026-05-26 12:45:50.899508: bru_input.bru_comp(z)
#> 2026-05-26 12:45:50.900805: bru_input.bm_pipe(z)
#> 2026-05-26 12:45:50.902279: bru_input.bm_multi(z:core)
#> 2026-05-26 12:45:50.903765: bru_input.bru_mapper(z:core:main)
#> 2026-05-26 12:45:50.905231: bru_input(bru_input) for (z)
#> 2026-05-26 12:45:50.912825: bru_input.bru_comp(Intercept)
#> 2026-05-26 12:45:50.914204: bru_input.bm_pipe(Intercept)
#> 2026-05-26 12:45:50.915679: bru_input.bm_multi(Intercept:core)
#> 2026-05-26 12:45:50.917124: bru_input.bru_mapper(Intercept:core:main)
#> 2026-05-26 12:45:50.918488: bru_input(bru_input) for (Intercept)
#> 2026-05-26 12:45:50.931396: iinla: Start
#> 2026-05-26 12:45:50.933124: iinla: Evaluate component linearisations
#> 2026-05-26 12:45:50.935338: Linearise components for each observation model
#> 2026-05-26 12:45:50.938335: Linearise component 'z'
#> 2026-05-26 12:45:50.944277: Linearise component 'Intercept'
#> 2026-05-26 12:45:50.947461: iinla: Evaluate component simplifications
#> 2026-05-26 12:45:50.948708: Simplify component mappers for each observation model
#> 2026-05-26 12:45:50.950175: Simplify component 'z'
#> 2026-05-26 12:45:50.953242: Simplify component 'Intercept'
#> 2026-05-26 12:45:50.956228: iinla: Evaluate predictor linearisation
#> 2026-05-26 12:45:50.95918: Linearise with respect to component 'z'
#> 2026-05-26 12:45:50.964411: Linearise with respect to component 'Intercept'
#> 2026-05-26 12:45:50.96786: iinla: Construct inla stack
#> 2026-05-26 12:45:50.9803: iinla: Model initialisation completed
#> 2026-05-26 12:45:50.981709: iinla: Iteration 1 [max: 10]
#> 2026-05-26 12:45:51.331385: iinla: Step rescaling: 61.8%, Contract (norm0 = 1969, norm1 = 1807, norm01 = 347.8)
#> 2026-05-26 12:45:51.333973: iinla: Step rescaling: 38.2%, Contract (norm0 = 398.9, norm1 = 295.1, norm01 = 347.8)
#> 2026-05-26 12:45:51.337423: iinla: Step rescaling: 27.4%, Approx Optimisation (norm0 = 184.1, norm1 = 225.5, norm01 = 347.8)
#> 2026-05-26 12:45:51.338928: iinla: |lin1-lin0| = 347.8
#>   <eta-lin1,delta>/|delta| = -198.3
#>   |eta-lin0 - delta <delta,eta-lin0>/<delta,delta>| = 107.4
#> 2026-05-26 12:45:51.340373: iinla: Step rescaling: 27.4% (norm0 = 184.1, norm1 = 225.5, norm01 = 347.8)
#> 2026-05-26 12:45:51.34205: iinla: Evaluate component linearisations
#> 2026-05-26 12:45:51.343146: Linearise components for each observation model
#> 2026-05-26 12:45:51.344445: Linearise component 'z'
#> 2026-05-26 12:45:51.347284: Linearise component 'Intercept'
#> 2026-05-26 12:45:51.349985: iinla: Evaluate predictor linearisation
#> 2026-05-26 12:45:51.35255: Linearise with respect to component 'z'
#> 2026-05-26 12:45:51.355917: Linearise with respect to component 'Intercept'
#> 2026-05-26 12:45:51.368339: iinla: Iteration 2 [max: 10]
#> 2026-05-26 12:45:51.94219: iinla: Step rescaling: 162%, Expand (norm0 = 365.1, norm1 = 141.6, norm01 = 225.4)
#> 2026-05-26 12:45:51.945187: iinla: Step rescaling: 100%, Overstep (norm0 = 225.4, norm1 = 10.71, norm01 = 225.4)
#> 2026-05-26 12:45:51.948639: iinla: Step rescaling: 99.69%, Approx Optimisation (norm0 = 224.7, norm1 = 10.68, norm01 = 225.4)
#> 2026-05-26 12:45:51.950029: iinla: |lin1-lin0| = 225.4
#>   <eta-lin1,delta>/|delta| = -1.008
#>   |eta-lin0 - delta <delta,eta-lin0>/<delta,delta>| = 10.63
#> 2026-05-26 12:45:51.95141: iinla: Step rescaling: 99.7% (norm0 = 224.7, norm1 = 10.68, norm01 = 225.4)
#> 2026-05-26 12:45:51.953083: iinla: Evaluate component linearisations
#> 2026-05-26 12:45:51.954207: Linearise components for each observation model
#> 2026-05-26 12:45:51.95552: Linearise component 'z'
#> 2026-05-26 12:45:51.958266: Linearise component 'Intercept'
#> 2026-05-26 12:45:51.961039: iinla: Evaluate predictor linearisation
#> 2026-05-26 12:45:51.963544: Linearise with respect to component 'z'
#> 2026-05-26 12:45:51.966935: Linearise with respect to component 'Intercept'
#> 2026-05-26 12:45:51.978858: iinla: Max deviation from previous: 46300% of SD, and line search is active
#> [stop if: < 10% and line search inactive]
#> 2026-05-26 12:45:51.981334: iinla: Iteration 3 [max: 10]
#> 2026-05-26 12:45:52.325184: iinla: Step rescaling: 162%, Expand (norm0 = 16.84, norm1 = 6.152, norm01 = 10.69)
#> 2026-05-26 12:45:52.327929: iinla: Step rescaling: 100%, Overstep (norm0 = 10.51, norm1 = 0.1742, norm01 = 10.69)
#> 2026-05-26 12:45:52.331352: iinla: Step rescaling: 101.7%, Approx Optimisation (norm0 = 10.69, norm1 = 0.01176, norm01 = 10.69)
#> 2026-05-26 12:45:52.332752: iinla: |lin1-lin0| = 10.69
#>   <eta-lin1,delta>/|delta| = -1.163e-05
#>   |eta-lin0 - delta <delta,eta-lin0>/<delta,delta>| = 0.01176
#> 2026-05-26 12:45:52.334209: iinla: Step rescaling: 102% (norm0 = 10.69, norm1 = 0.01176, norm01 = 10.69)
#> 2026-05-26 12:45:52.335913: iinla: Evaluate component linearisations
#> 2026-05-26 12:45:52.33718: Linearise components for each observation model
#> 2026-05-26 12:45:52.338687: Linearise component 'z'
#> 2026-05-26 12:45:52.341873: Linearise component 'Intercept'
#> 2026-05-26 12:45:52.344736: iinla: Evaluate predictor linearisation
#> 2026-05-26 12:45:52.347479: Linearise with respect to component 'z'
#> 2026-05-26 12:45:52.351111: Linearise with respect to component 'Intercept'
#> 2026-05-26 12:45:52.363075: iinla: Max deviation from previous: 496% of SD, and line search is active
#> [stop if: < 10% and line search inactive]
#> 2026-05-26 12:45:52.365527: iinla: Iteration 4 [max: 10]
#> 2026-05-26 12:45:52.734427: iinla: Step rescaling: 162%, Expand (norm0 = 0.01903, norm1 = 0.007269, norm01 = 0.01176)
#> 2026-05-26 12:45:52.737399: iinla: Step rescaling: 100%, Overstep (norm0 = 0.01176, norm1 = 5.869e-08, norm01 = 0.01176)
#> 2026-05-26 12:45:52.742109: iinla: Step rescaling: 100%, Approx Optimisation (norm0 = 0.01176, norm1 = 5.856e-08, norm01 = 0.01176)
#> 2026-05-26 12:45:52.743824: iinla: |lin1-lin0| = 0.01176
#>   <eta-lin1,delta>/|delta| = 8.233e-11
#>   |eta-lin0 - delta <delta,eta-lin0>/<delta,delta>| = 5.856e-08
#> 2026-05-26 12:45:52.745797: iinla: Evaluate component linearisations
#> 2026-05-26 12:45:52.747015: Linearise components for each observation model
#> 2026-05-26 12:45:52.748457: Linearise component 'z'
#> 2026-05-26 12:45:52.751584: Linearise component 'Intercept'
#> 2026-05-26 12:45:52.754388: iinla: Evaluate predictor linearisation
#> 2026-05-26 12:45:52.757036: Linearise with respect to component 'z'
#> 2026-05-26 12:45:52.760485: Linearise with respect to component 'Intercept'
#> 2026-05-26 12:45:52.773334: iinla: Max deviation from previous: 8.06% of SD, and line search is inactive
#> [stop if: < 10% and line search inactive]
#> 2026-05-26 12:45:52.77474: iinla: Convergence criterion met.
#>        Running final INLA integration step with known theta mode.
#> 2026-05-26 12:45:52.777451: iinla: Iteration 5 [max: 10]
#> 2026-05-26 12:45:53.137992: iinla: Computation completed
#> 2026-05-26 12:45:54.292735: bru: Preprocessing
#> 2026-05-26 12:45:54.302197: Evaluate component inputs for each observation model
#> 2026-05-26 12:45:54.303517: bru_input(bru_comp_list)
#> 2026-05-26 12:45:54.305092: bru_input.bru_comp(field)
#> 2026-05-26 12:45:54.30652: bru_input.bm_pipe(field)
#> 2026-05-26 12:45:54.307995: bru_input.bm_multi(field:core)
#> 2026-05-26 12:45:54.309449: bru_input.bru_mapper(field:core:main)
#> 2026-05-26 12:45:54.310755: bru_input(bru_input) for (field)
#> 2026-05-26 12:45:54.323472: iinla: Start
#> 2026-05-26 12:45:54.325071: iinla: Evaluate component linearisations
#> 2026-05-26 12:45:54.326224: Linearise components for each observation model
#> 2026-05-26 12:45:54.327564: Linearise component 'field'
#> 2026-05-26 12:45:54.334477: iinla: Evaluate component simplifications
#> 2026-05-26 12:45:54.335758: Simplify component mappers for each observation model
#> 2026-05-26 12:45:54.337274: Simplify component 'field'
#> 2026-05-26 12:45:54.344011: iinla: Evaluate predictor linearisation
#> 2026-05-26 12:45:54.3474: iinla: Construct inla stack
#> 2026-05-26 12:45:54.357513: iinla: Model initialisation completed
#> 2026-05-26 12:45:54.358945: iinla: Iteration 1 [max: 1]
#> 2026-05-26 12:45:54.837335: iinla: Computation completed
#> 2026-05-26 12:45:56.85779: bru_input(bru_comp_list)
#> 2026-05-26 12:45:59.674772: bru: Preprocessing
#> 2026-05-26 12:45:59.68834: Evaluate component inputs for each observation model
#> 2026-05-26 12:45:59.689585: bru_input(bru_comp_list)
#> 2026-05-26 12:45:59.690995: bru_input.bru_comp(x)
#> 2026-05-26 12:45:59.692243: bru_input.bm_pipe(x)
#> 2026-05-26 12:45:59.693641: bru_input.bm_multi(x:core)
#> 2026-05-26 12:45:59.695051: bru_input.bru_mapper(x:core:main)
#> 2026-05-26 12:45:59.69638: bru_input(bru_input) for (x)
#> 2026-05-26 12:45:59.703833: bru_input(bru_comp_list)
#> 2026-05-26 12:45:59.705882: bru_input.bru_comp(x)
#> 2026-05-26 12:45:59.70743: bru_input.bm_pipe(x)
#> 2026-05-26 12:45:59.714772: bru_input.bm_multi(x:core)
#> 2026-05-26 12:45:59.716385: bru_input.bru_mapper(x:core:main)
#> 2026-05-26 12:45:59.717792: bru_input(bru_input) for (x)
#> 2026-05-26 12:46:01.82925: bru_input(bru_input) for (LABEL)
#> 2026-05-26 12:46:03.713366: bru_input_text(bru_input) for (LABEL)
#> 2026-05-26 12:46:03.723452: bru_input.bru_comp(x)
#> 2026-05-26 12:46:03.724605: bru_input_text.bm_pipe(x)
#> 2026-05-26 12:46:03.725761: bru_input_text.bm_multi(x:core)
#> 2026-05-26 12:46:03.726888: bru_input_text.bru_mapper(x:core:main)
#> 2026-05-26 12:46:03.727983: bru_input_text(bru_input) for (x)

bru_log(verbosity = 2L)
#> 2026-05-26 12:45:38.190956: inlabru loaded
#> 2026-05-26 12:45:38.191377: Clear override options
#> 2026-05-26 12:45:49.727601: bru: Preprocessing
#> 2026-05-26 12:45:49.84122: iinla: Iteration 1 [max: 1]
#> 2026-05-26 12:45:50.418018: bru: Preprocessing
#> 2026-05-26 12:45:50.504482: iinla: Iteration 1 [max: 1]
#> 2026-05-26 12:45:50.888623: bru: Preprocessing
#> 2026-05-26 12:45:50.981709: iinla: Iteration 1 [max: 10]
#> 2026-05-26 12:45:51.340373: iinla: Step rescaling: 27.4% (norm0 = 184.1, norm1 = 225.5, norm01 = 347.8)
#> 2026-05-26 12:45:51.368339: iinla: Iteration 2 [max: 10]
#> 2026-05-26 12:45:51.95141: iinla: Step rescaling: 99.7% (norm0 = 224.7, norm1 = 10.68, norm01 = 225.4)
#> 2026-05-26 12:45:51.978858: iinla: Max deviation from previous: 46300% of SD, and line search is active
#> [stop if: < 10% and line search inactive]
#> 2026-05-26 12:45:51.981334: iinla: Iteration 3 [max: 10]
#> 2026-05-26 12:45:52.334209: iinla: Step rescaling: 102% (norm0 = 10.69, norm1 = 0.01176, norm01 = 10.69)
#> 2026-05-26 12:45:52.363075: iinla: Max deviation from previous: 496% of SD, and line search is active
#> [stop if: < 10% and line search inactive]
#> 2026-05-26 12:45:52.365527: iinla: Iteration 4 [max: 10]
#> 2026-05-26 12:45:52.773334: iinla: Max deviation from previous: 8.06% of SD, and line search is inactive
#> [stop if: < 10% and line search inactive]
#> 2026-05-26 12:45:52.77474: iinla: Convergence criterion met.
#>        Running final INLA integration step with known theta mode.
#> 2026-05-26 12:45:52.777451: iinla: Iteration 5 [max: 10]
#> 2026-05-26 12:45:54.292735: bru: Preprocessing
#> 2026-05-26 12:45:54.358945: iinla: Iteration 1 [max: 1]
#> 2026-05-26 12:45:59.674772: bru: Preprocessing
print(bru_log(), timestamp = TRUE, verbosity = TRUE)
#> 2026-05-26 12:45:38.190956: inlabru loaded (level 1)
#> 2026-05-26 12:45:38.191377: Clear override options (level 1)
#> 2026-05-26 12:45:49.727601: bru: Preprocessing (level 1)
#> 2026-05-26 12:45:49.743866: Evaluate component inputs for each observation model (level 3)
#> 2026-05-26 12:45:49.745247: bru_input(bru_comp_list) (level 4)
#> 2026-05-26 12:45:49.746636: bru_input.bru_comp(x) (level 4)
#> 2026-05-26 12:45:49.748101: bru_input.bm_pipe(x) (level 5)
#> 2026-05-26 12:45:49.749836: bru_input.bm_multi(x:core) (level 5)
#> 2026-05-26 12:45:49.751425: bru_input.bru_mapper(x:core:main) (level 5)
#> 2026-05-26 12:45:49.753049: bru_input(bru_input) for (x) (level 5)
#> 2026-05-26 12:45:49.761505: bru_input.bru_comp(Intercept) (level 4)
#> 2026-05-26 12:45:49.762754: bru_input.bm_pipe(Intercept) (level 5)
#> 2026-05-26 12:45:49.764135: bru_input.bm_multi(Intercept:core) (level 5)
#> 2026-05-26 12:45:49.765461: bru_input.bru_mapper(Intercept:core:main) (level 5)
#> 2026-05-26 12:45:49.766691: bru_input(bru_input) for (Intercept) (level 5)
#> 2026-05-26 12:45:49.78371: iinla: Start (level 3)
#> 2026-05-26 12:45:49.785265: iinla: Evaluate component linearisations (level 3)
#> 2026-05-26 12:45:49.786408: Linearise components for each observation model (level 3)
#> 2026-05-26 12:45:49.787807: Linearise component 'x' (level 4)
#> 2026-05-26 12:45:49.790721: Linearise component 'Intercept' (level 4)
#> 2026-05-26 12:45:49.793379: iinla: Evaluate component simplifications (level 3)
#> 2026-05-26 12:45:49.794529: Simplify component mappers for each observation model (level 3)
#> 2026-05-26 12:45:49.795882: Simplify component 'x' (level 4)
#> 2026-05-26 12:45:49.798847: Simplify component 'Intercept' (level 4)
#> 2026-05-26 12:45:49.801742: iinla: Evaluate predictor linearisation (level 3)
#> 2026-05-26 12:45:49.815638: iinla: Construct inla stack (level 3)
#> 2026-05-26 12:45:49.839587: iinla: Model initialisation completed (level 3)
#> 2026-05-26 12:45:49.84122: iinla: Iteration 1 [max: 1] (level 1)
#> 2026-05-26 12:45:50.39276: iinla: Computation completed (level 3)
#> 2026-05-26 12:45:50.418018: bru: Preprocessing (level 1)
#> 2026-05-26 12:45:50.426873: Evaluate component inputs for each observation model (level 3)
#> 2026-05-26 12:45:50.428097: bru_input(bru_comp_list) (level 4)
#> 2026-05-26 12:45:50.429435: bru_input.bru_comp(x) (level 4)
#> 2026-05-26 12:45:50.430851: bru_input.bm_pipe(x) (level 5)
#> 2026-05-26 12:45:50.432316: bru_input.bm_multi(x:core) (level 5)
#> 2026-05-26 12:45:50.433751: bru_input.bru_mapper(x:core:main) (level 5)
#> 2026-05-26 12:45:50.435068: bru_input(bru_input) for (x) (level 5)
#> 2026-05-26 12:45:50.445308: bru_input.bru_comp(Intercept) (level 4)
#> 2026-05-26 12:45:50.446794: bru_input.bm_pipe(Intercept) (level 5)
#> 2026-05-26 12:45:50.448344: bru_input.bm_multi(Intercept:core) (level 5)
#> 2026-05-26 12:45:50.44979: bru_input.bru_mapper(Intercept:core:main) (level 5)
#> 2026-05-26 12:45:50.451172: bru_input(bru_input) for (Intercept) (level 5)
#> 2026-05-26 12:45:50.464854: iinla: Start (level 3)
#> 2026-05-26 12:45:50.466457: iinla: Evaluate component linearisations (level 3)
#> 2026-05-26 12:45:50.467615: Linearise components for each observation model (level 3)
#> 2026-05-26 12:45:50.46896: Linearise component 'x' (level 4)
#> 2026-05-26 12:45:50.472423: Linearise component 'Intercept' (level 4)
#> 2026-05-26 12:45:50.475768: iinla: Evaluate component simplifications (level 3)
#> 2026-05-26 12:45:50.477107: Simplify component mappers for each observation model (level 3)
#> 2026-05-26 12:45:50.47861: Simplify component 'x' (level 4)
#> 2026-05-26 12:45:50.481707: Simplify component 'Intercept' (level 4)
#> 2026-05-26 12:45:50.484711: iinla: Evaluate predictor linearisation (level 3)
#> 2026-05-26 12:45:50.491022: iinla: Construct inla stack (level 3)
#> 2026-05-26 12:45:50.502785: iinla: Model initialisation completed (level 3)
#> 2026-05-26 12:45:50.504482: iinla: Iteration 1 [max: 1] (level 1)
#> 2026-05-26 12:45:50.867193: iinla: Computation completed (level 3)
#> 2026-05-26 12:45:50.888623: bru: Preprocessing (level 1)
#> 2026-05-26 12:45:50.896928: Evaluate component inputs for each observation model (level 3)
#> 2026-05-26 12:45:50.898171: bru_input(bru_comp_list) (level 4)
#> 2026-05-26 12:45:50.899508: bru_input.bru_comp(z) (level 4)
#> 2026-05-26 12:45:50.900805: bru_input.bm_pipe(z) (level 5)
#> 2026-05-26 12:45:50.902279: bru_input.bm_multi(z:core) (level 5)
#> 2026-05-26 12:45:50.903765: bru_input.bru_mapper(z:core:main) (level 5)
#> 2026-05-26 12:45:50.905231: bru_input(bru_input) for (z) (level 5)
#> 2026-05-26 12:45:50.912825: bru_input.bru_comp(Intercept) (level 4)
#> 2026-05-26 12:45:50.914204: bru_input.bm_pipe(Intercept) (level 5)
#> 2026-05-26 12:45:50.915679: bru_input.bm_multi(Intercept:core) (level 5)
#> 2026-05-26 12:45:50.917124: bru_input.bru_mapper(Intercept:core:main) (level 5)
#> 2026-05-26 12:45:50.918488: bru_input(bru_input) for (Intercept) (level 5)
#> 2026-05-26 12:45:50.931396: iinla: Start (level 3)
#> 2026-05-26 12:45:50.933124: iinla: Evaluate component linearisations (level 3)
#> 2026-05-26 12:45:50.935338: Linearise components for each observation model (level 3)
#> 2026-05-26 12:45:50.938335: Linearise component 'z' (level 4)
#> 2026-05-26 12:45:50.944277: Linearise component 'Intercept' (level 4)
#> 2026-05-26 12:45:50.947461: iinla: Evaluate component simplifications (level 3)
#> 2026-05-26 12:45:50.948708: Simplify component mappers for each observation model (level 3)
#> 2026-05-26 12:45:50.950175: Simplify component 'z' (level 4)
#> 2026-05-26 12:45:50.953242: Simplify component 'Intercept' (level 4)
#> 2026-05-26 12:45:50.956228: iinla: Evaluate predictor linearisation (level 3)
#> 2026-05-26 12:45:50.95918: Linearise with respect to component 'z' (level 5)
#> 2026-05-26 12:45:50.964411: Linearise with respect to component 'Intercept' (level 5)
#> 2026-05-26 12:45:50.96786: iinla: Construct inla stack (level 3)
#> 2026-05-26 12:45:50.9803: iinla: Model initialisation completed (level 3)
#> 2026-05-26 12:45:50.981709: iinla: Iteration 1 [max: 10] (level 1)
#> 2026-05-26 12:45:51.331385: iinla: Step rescaling: 61.8%, Contract (norm0 = 1969, norm1 = 1807, norm01 = 347.8) (level 3)
#> 2026-05-26 12:45:51.333973: iinla: Step rescaling: 38.2%, Contract (norm0 = 398.9, norm1 = 295.1, norm01 = 347.8) (level 3)
#> 2026-05-26 12:45:51.337423: iinla: Step rescaling: 27.4%, Approx Optimisation (norm0 = 184.1, norm1 = 225.5, norm01 = 347.8) (level 3)
#> 2026-05-26 12:45:51.338928: iinla: |lin1-lin0| = 347.8
#>   <eta-lin1,delta>/|delta| = -198.3
#>   |eta-lin0 - delta <delta,eta-lin0>/<delta,delta>| = 107.4 (level 4)
#> 2026-05-26 12:45:51.340373: iinla: Step rescaling: 27.4% (norm0 = 184.1, norm1 = 225.5, norm01 = 347.8) (level 2)
#> 2026-05-26 12:45:51.34205: iinla: Evaluate component linearisations (level 3)
#> 2026-05-26 12:45:51.343146: Linearise components for each observation model (level 3)
#> 2026-05-26 12:45:51.344445: Linearise component 'z' (level 4)
#> 2026-05-26 12:45:51.347284: Linearise component 'Intercept' (level 4)
#> 2026-05-26 12:45:51.349985: iinla: Evaluate predictor linearisation (level 3)
#> 2026-05-26 12:45:51.35255: Linearise with respect to component 'z' (level 5)
#> 2026-05-26 12:45:51.355917: Linearise with respect to component 'Intercept' (level 5)
#> 2026-05-26 12:45:51.368339: iinla: Iteration 2 [max: 10] (level 1)
#> 2026-05-26 12:45:51.94219: iinla: Step rescaling: 162%, Expand (norm0 = 365.1, norm1 = 141.6, norm01 = 225.4) (level 3)
#> 2026-05-26 12:45:51.945187: iinla: Step rescaling: 100%, Overstep (norm0 = 225.4, norm1 = 10.71, norm01 = 225.4) (level 3)
#> 2026-05-26 12:45:51.948639: iinla: Step rescaling: 99.69%, Approx Optimisation (norm0 = 224.7, norm1 = 10.68, norm01 = 225.4) (level 3)
#> 2026-05-26 12:45:51.950029: iinla: |lin1-lin0| = 225.4
#>   <eta-lin1,delta>/|delta| = -1.008
#>   |eta-lin0 - delta <delta,eta-lin0>/<delta,delta>| = 10.63 (level 4)
#> 2026-05-26 12:45:51.95141: iinla: Step rescaling: 99.7% (norm0 = 224.7, norm1 = 10.68, norm01 = 225.4) (level 2)
#> 2026-05-26 12:45:51.953083: iinla: Evaluate component linearisations (level 3)
#> 2026-05-26 12:45:51.954207: Linearise components for each observation model (level 3)
#> 2026-05-26 12:45:51.95552: Linearise component 'z' (level 4)
#> 2026-05-26 12:45:51.958266: Linearise component 'Intercept' (level 4)
#> 2026-05-26 12:45:51.961039: iinla: Evaluate predictor linearisation (level 3)
#> 2026-05-26 12:45:51.963544: Linearise with respect to component 'z' (level 5)
#> 2026-05-26 12:45:51.966935: Linearise with respect to component 'Intercept' (level 5)
#> 2026-05-26 12:45:51.978858: iinla: Max deviation from previous: 46300% of SD, and line search is active
#> [stop if: < 10% and line search inactive] (level 1)
#> 2026-05-26 12:45:51.981334: iinla: Iteration 3 [max: 10] (level 1)
#> 2026-05-26 12:45:52.325184: iinla: Step rescaling: 162%, Expand (norm0 = 16.84, norm1 = 6.152, norm01 = 10.69) (level 3)
#> 2026-05-26 12:45:52.327929: iinla: Step rescaling: 100%, Overstep (norm0 = 10.51, norm1 = 0.1742, norm01 = 10.69) (level 3)
#> 2026-05-26 12:45:52.331352: iinla: Step rescaling: 101.7%, Approx Optimisation (norm0 = 10.69, norm1 = 0.01176, norm01 = 10.69) (level 3)
#> 2026-05-26 12:45:52.332752: iinla: |lin1-lin0| = 10.69
#>   <eta-lin1,delta>/|delta| = -1.163e-05
#>   |eta-lin0 - delta <delta,eta-lin0>/<delta,delta>| = 0.01176 (level 4)
#> 2026-05-26 12:45:52.334209: iinla: Step rescaling: 102% (norm0 = 10.69, norm1 = 0.01176, norm01 = 10.69) (level 2)
#> 2026-05-26 12:45:52.335913: iinla: Evaluate component linearisations (level 3)
#> 2026-05-26 12:45:52.33718: Linearise components for each observation model (level 3)
#> 2026-05-26 12:45:52.338687: Linearise component 'z' (level 4)
#> 2026-05-26 12:45:52.341873: Linearise component 'Intercept' (level 4)
#> 2026-05-26 12:45:52.344736: iinla: Evaluate predictor linearisation (level 3)
#> 2026-05-26 12:45:52.347479: Linearise with respect to component 'z' (level 5)
#> 2026-05-26 12:45:52.351111: Linearise with respect to component 'Intercept' (level 5)
#> 2026-05-26 12:45:52.363075: iinla: Max deviation from previous: 496% of SD, and line search is active
#> [stop if: < 10% and line search inactive] (level 1)
#> 2026-05-26 12:45:52.365527: iinla: Iteration 4 [max: 10] (level 1)
#> 2026-05-26 12:45:52.734427: iinla: Step rescaling: 162%, Expand (norm0 = 0.01903, norm1 = 0.007269, norm01 = 0.01176) (level 3)
#> 2026-05-26 12:45:52.737399: iinla: Step rescaling: 100%, Overstep (norm0 = 0.01176, norm1 = 5.869e-08, norm01 = 0.01176) (level 3)
#> 2026-05-26 12:45:52.742109: iinla: Step rescaling: 100%, Approx Optimisation (norm0 = 0.01176, norm1 = 5.856e-08, norm01 = 0.01176) (level 3)
#> 2026-05-26 12:45:52.743824: iinla: |lin1-lin0| = 0.01176
#>   <eta-lin1,delta>/|delta| = 8.233e-11
#>   |eta-lin0 - delta <delta,eta-lin0>/<delta,delta>| = 5.856e-08 (level 4)
#> 2026-05-26 12:45:52.745797: iinla: Evaluate component linearisations (level 3)
#> 2026-05-26 12:45:52.747015: Linearise components for each observation model (level 3)
#> 2026-05-26 12:45:52.748457: Linearise component 'z' (level 4)
#> 2026-05-26 12:45:52.751584: Linearise component 'Intercept' (level 4)
#> 2026-05-26 12:45:52.754388: iinla: Evaluate predictor linearisation (level 3)
#> 2026-05-26 12:45:52.757036: Linearise with respect to component 'z' (level 5)
#> 2026-05-26 12:45:52.760485: Linearise with respect to component 'Intercept' (level 5)
#> 2026-05-26 12:45:52.773334: iinla: Max deviation from previous: 8.06% of SD, and line search is inactive
#> [stop if: < 10% and line search inactive] (level 1)
#> 2026-05-26 12:45:52.77474: iinla: Convergence criterion met.
#>        Running final INLA integration step with known theta mode. (level 1)
#> 2026-05-26 12:45:52.777451: iinla: Iteration 5 [max: 10] (level 1)
#> 2026-05-26 12:45:53.137992: iinla: Computation completed (level 3)
#> 2026-05-26 12:45:54.292735: bru: Preprocessing (level 1)
#> 2026-05-26 12:45:54.302197: Evaluate component inputs for each observation model (level 3)
#> 2026-05-26 12:45:54.303517: bru_input(bru_comp_list) (level 4)
#> 2026-05-26 12:45:54.305092: bru_input.bru_comp(field) (level 4)
#> 2026-05-26 12:45:54.30652: bru_input.bm_pipe(field) (level 5)
#> 2026-05-26 12:45:54.307995: bru_input.bm_multi(field:core) (level 5)
#> 2026-05-26 12:45:54.309449: bru_input.bru_mapper(field:core:main) (level 5)
#> 2026-05-26 12:45:54.310755: bru_input(bru_input) for (field) (level 5)
#> 2026-05-26 12:45:54.323472: iinla: Start (level 3)
#> 2026-05-26 12:45:54.325071: iinla: Evaluate component linearisations (level 3)
#> 2026-05-26 12:45:54.326224: Linearise components for each observation model (level 3)
#> 2026-05-26 12:45:54.327564: Linearise component 'field' (level 4)
#> 2026-05-26 12:45:54.334477: iinla: Evaluate component simplifications (level 3)
#> 2026-05-26 12:45:54.335758: Simplify component mappers for each observation model (level 3)
#> 2026-05-26 12:45:54.337274: Simplify component 'field' (level 4)
#> 2026-05-26 12:45:54.344011: iinla: Evaluate predictor linearisation (level 3)
#> 2026-05-26 12:45:54.3474: iinla: Construct inla stack (level 3)
#> 2026-05-26 12:45:54.357513: iinla: Model initialisation completed (level 3)
#> 2026-05-26 12:45:54.358945: iinla: Iteration 1 [max: 1] (level 1)
#> 2026-05-26 12:45:54.837335: iinla: Computation completed (level 3)
#> 2026-05-26 12:45:56.85779: bru_input(bru_comp_list) (level 4)
#> 2026-05-26 12:45:59.674772: bru: Preprocessing (level 1)
#> 2026-05-26 12:45:59.68834: Evaluate component inputs for each observation model (level 3)
#> 2026-05-26 12:45:59.689585: bru_input(bru_comp_list) (level 4)
#> 2026-05-26 12:45:59.690995: bru_input.bru_comp(x) (level 4)
#> 2026-05-26 12:45:59.692243: bru_input.bm_pipe(x) (level 5)
#> 2026-05-26 12:45:59.693641: bru_input.bm_multi(x:core) (level 5)
#> 2026-05-26 12:45:59.695051: bru_input.bru_mapper(x:core:main) (level 5)
#> 2026-05-26 12:45:59.69638: bru_input(bru_input) for (x) (level 5)
#> 2026-05-26 12:45:59.703833: bru_input(bru_comp_list) (level 4)
#> 2026-05-26 12:45:59.705882: bru_input.bru_comp(x) (level 4)
#> 2026-05-26 12:45:59.70743: bru_input.bm_pipe(x) (level 5)
#> 2026-05-26 12:45:59.714772: bru_input.bm_multi(x:core) (level 5)
#> 2026-05-26 12:45:59.716385: bru_input.bru_mapper(x:core:main) (level 5)
#> 2026-05-26 12:45:59.717792: bru_input(bru_input) for (x) (level 5)
#> 2026-05-26 12:46:01.82925: bru_input(bru_input) for (LABEL) (level 5)
#> 2026-05-26 12:46:03.713366: bru_input_text(bru_input) for (LABEL) (level 5)
#> 2026-05-26 12:46:03.723452: bru_input.bru_comp(x) (level 4)
#> 2026-05-26 12:46:03.724605: bru_input_text.bm_pipe(x) (level 5)
#> 2026-05-26 12:46:03.725761: bru_input_text.bm_multi(x:core) (level 5)
#> 2026-05-26 12:46:03.726888: bru_input_text.bru_mapper(x:core:main) (level 5)
#> 2026-05-26 12:46:03.727983: bru_input_text(bru_input) for (x) (level 5)
```
