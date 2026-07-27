# Additional bru options

Construct a `bru_options` object including the default and global
options, and converting deprecated option names.

## Usage

``` r
bru_call_options(...)
```

## Arguments

- ...:

  Options passed on to
  [`as.bru_options()`](https://inlabru-org.github.io/inlabru/reference/bru_options.md)

## Value

A `bru_options` object

## Author

Finn Lindgren <finn.lindgren@gmail.com>

## Examples

``` r
# \donttest{

opts <- bru_call_options()

# Print them:
opts
#> $bru_verbose
#> [1] 0
#> 
#> $bru_verbose_store
#> [1] Inf
#> 
#> $bru_max_iter
#> [1] 10
#> 
#> $bru_run
#> [1] TRUE
#> 
#> $bru_int_args
#> $bru_int_args$method
#> [1] "stable"
#> 
#> $bru_int_args$nsub1
#> [1] 30
#> 
#> $bru_int_args$nsub2
#> [1] 9
#> 
#> 
#> $bru_method
#> $bru_method$autodiff
#> [1] "fullchain"
#> 
#> $bru_method$agg
#> [1] "fullchain"
#> 
#> $bru_method$finite_diff
#> [1] "forward"
#> 
#> $bru_method$taylor
#> [1] "pandemic"
#> 
#> $bru_method$search
#> [1] "all"
#> 
#> $bru_method$factor
#> [1] 1.618034
#> 
#> $bru_method$rel_tol
#> [1] 0.1
#> 
#> $bru_method$max_step
#> [1] 2
#> 
#> $bru_method$line_opt_method
#> [1] "onestep"
#> 
#> 
#> $bru_compress_cp
#> [1] FALSE
#> 
#> $bru_debug
#> [1] FALSE
#> 
#> $bru_compat_pre_2_14_enable
#> [1] TRUE
#> 
#> $E
#> [1] 1
#> 
#> $Ntrials
#> [1] 1
#> 
#> $control.compute
#> $control.compute$config
#> [1] TRUE
#> 
#> $control.compute$control.gcpo
#> list()
#> 
#> 
#> $control.inla
#> $control.inla$int.strategy
#> [1] "auto"
#> 
#> 
#> $control.fixed
#> $control.fixed$expand.factor.strategy
#> [1] "inla"
#> 
#> 
#> attr(,"class")
#> [1] "bru_options" "list"       
# }
```
