# Mapper for general expressions

Constructs an expression mapping

## Usage

``` r
bm_expr(
  expr,
  assume = character(0),
  labels = NULL,
  .envir = rlang::caller_env()
)
```

## Arguments

- expr:

  An expression object, typically created using
  [`rlang::quo()`](https://rlang.r-lib.org/reference/defusing-advanced.html),
  that can be evaluated in a data mask containing the root variables,
  derived variables, and data variables.

- assume:

  character vector listing valid assumptions for the combination of the
  expression and data, derived variables, and root variables. This can
  be used to specify assumptions that allow for more efficient Jacobian
  calculations. For example, if the expression is row-wise in the data
  and derived variables, then `assume = "rowwise"` can be used to
  indicate that the Jacobian with respect to derived variables can be
  calculated simultaneously for all rows. Supported assumptions are
  "rowwise", "linear", "additive", and "no_root".

- labels:

  list of `root`, `derived`, and `suffix`. The `root` and `derived`
  elements specify the pronoun labels for the data mask for root and
  derived variables. Default "root" and "derived", respectively. For
  example, if `root = "latent"`, then the root variables will be
  accessible in the expression as `.latent$varname`. If `suffix` is
  non-NULL, a character string specifying a suffix to add to the root
  variable names when making them directly available to the expression.
  Default is NULL, equivalent to `""`. For example, if
  `suffix = "_latent"`, then a root variable named `"x"` will be
  available as `"x_latent"` in the expression, in addition to the data
  mask pronoun version.

- .envir:

  The environment for the expression evaluation. By default, this is set
  to the caller environment.

## Value

A `bm_expr` mapper object

## Details

The `input` is currently ignored.

`data` should be a list with data objects, with the main object called
`data`. If `is_rowwise == TRUE`, the number of rows in the `data`
data.frame determines the number of rows in the output, and the columns
can be used as constants in the expression, accessed via
`.data$colname`, `.data.$colname`, or `.data.[["colname"]]`.

## See also

[bru_mapper](https://inlabru-org.github.io/inlabru/reference/bru_mapper.md),
[bru_mapper_generics](https://inlabru-org.github.io/inlabru/reference/bru_mapper_generics.md)

Other mappers:
[`bm_aggregate()`](https://inlabru-org.github.io/inlabru/reference/bm_aggregate.md),
[`bm_collect()`](https://inlabru-org.github.io/inlabru/reference/bm_collect.md),
[`bm_const()`](https://inlabru-org.github.io/inlabru/reference/bm_const.md),
[`bm_factor()`](https://inlabru-org.github.io/inlabru/reference/bm_factor.md),
[`bm_fm_mesh_1d`](https://inlabru-org.github.io/inlabru/reference/bm_fm_mesh_1d.md),
[`bm_fmesher()`](https://inlabru-org.github.io/inlabru/reference/bm_fmesher.md),
[`bm_harmonics()`](https://inlabru-org.github.io/inlabru/reference/bm_harmonics.md),
[`bm_index()`](https://inlabru-org.github.io/inlabru/reference/bm_index.md),
[`bm_linear()`](https://inlabru-org.github.io/inlabru/reference/bm_linear.md),
[`bm_logitaverage()`](https://inlabru-org.github.io/inlabru/reference/bm_logitaverage.md),
[`bm_logsumexp()`](https://inlabru-org.github.io/inlabru/reference/bm_logsumexp.md),
[`bm_marginal()`](https://inlabru-org.github.io/inlabru/reference/bm_marginal.md),
[`bm_matrix()`](https://inlabru-org.github.io/inlabru/reference/bm_matrix.md),
[`bm_multi()`](https://inlabru-org.github.io/inlabru/reference/bm_multi.md),
[`bm_pipe()`](https://inlabru-org.github.io/inlabru/reference/bm_pipe.md),
[`bm_reparam()`](https://inlabru-org.github.io/inlabru/reference/bm_reparam.md),
[`bm_repeat()`](https://inlabru-org.github.io/inlabru/reference/bm_repeat.md),
[`bm_scale()`](https://inlabru-org.github.io/inlabru/reference/bm_scale.md),
[`bm_shift()`](https://inlabru-org.github.io/inlabru/reference/bm_shift.md),
[`bm_sum()`](https://inlabru-org.github.io/inlabru/reference/bm_sum.md),
[`bm_taylor()`](https://inlabru-org.github.io/inlabru/reference/bm_taylor.md),
[`bru_get_mapper()`](https://inlabru-org.github.io/inlabru/reference/bru_get_mapper.md),
[`bru_mapper()`](https://inlabru-org.github.io/inlabru/reference/bru_mapper.md)

## Examples

``` r
# Basic expression with only root variables ("x").
(m <- bm_expr(rlang::quo(cos(x))))
#> expr(~cos(x))
ibm_eval(m, list(), list(x = 1:5))
#> [1]  0.5403023 -0.4161468 -0.9899925 -0.6536436  0.2836622
ibm_eval2(m, list(), list(x = 1:5))
#> $offset
#> [1]  0.5403023 -0.4161468 -0.9899925 -0.6536436  0.2836622
#> 
#> $jacobian
#> 5 x 5 sparse Matrix of class "dgCMatrix"
#>                                                          
#> [1,] -0.8414713  .          .         .         .        
#> [2,]  .         -0.9092972  .         .         .        
#> [3,]  .          .         -0.1411195 .         .        
#> [4,]  .          .          .         0.7568028 .        
#> [5,]  .          .          .         .         0.9589241
#> 

# Expression with data
(m <- bm_expr(rlang::quo(cos(x) * .data$z)))
#> expr(~cos(x) * .data$z)
the_data <- list(data = data.frame(z = 11:15))
ibm_eval(m, list(), list(x = 1:5), data = the_data)
#> [1]   5.943325  -4.993762 -12.869902  -9.151011   4.254933
ibm_eval2(m, list(), list(x = 1:5), data = the_data)
#> $offset
#> [1]   5.943325  -4.993762 -12.869902  -9.151011   4.254933
#> 
#> $jacobian
#> 5 x 5 sparse Matrix of class "dgCMatrix"
#>                                                     
#> [1,] -9.256184   .        .         .        .      
#> [2,]  .        -10.91157  .         .        .      
#> [3,]  .          .       -1.834554  .        .      
#> [4,]  .          .        .        10.59524  .      
#> [5,]  .          .        .         .       14.38386
#> 

# Expression with data, root variables, and derived variables.
(m <- bm_expr(
  rlang::quo(sin(x_latent) + cos(.effect$y) * .data$z),
  labels = list(root = "latent", derived = "effect", suffix = "_latent")
)
)
#> expr(~sin(x_latent) + cos(.effect$y) * .data$z)
ibm_eval(
  m,
  list(),
  derived = list(y = 2:6), # y = x + 1
  jacobians = list(y = list(x = Matrix::Diagonal(1.0, 5))),
  data = the_data,
  state = list(x = 1:5)
)
#> [1]  -3.736144 -10.970613  -8.356247   3.214468  13.443630
```
