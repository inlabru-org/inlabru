# Methods for mapper lists

`bru_mapper` lists can be combined into `bm_list` lists.

## Usage

``` r
as_bm_list(x)

# S3 method for class 'list'
as_bm_list(x)

# S3 method for class 'bm_list'
as_bm_list(x)

# S3 method for class 'bru_comp_list'
as_bm_list(x)

# S3 method for class 'bru_mapper'
c(...)

# S3 method for class 'bm_list'
c(...)

# S3 method for class 'bm_list'
x[i]
```

## Arguments

- x:

  `bm_list` object from which to extract element(s)

- ...:

  Objects to be combined.

- i:

  indices specifying elements to extract

## Value

A `bm_list` object

## Methods (by generic)

- `c(bm_list)`: The `...` arguments should be `bm_list` objects.

- `[`: Extract sub-list

## Functions

- `c(bru_mapper)`: The `...` arguments should be `bru_mapper` objects.

## Examples

``` r
m <- c(A = bm_const(), B = bm_scale())
str(m)
#> List of 2
#>  $ A: list()
#>   ..- attr(*, "class")= chr [1:2] "bm_const" "bru_mapper"
#>  $ B: list()
#>   ..- attr(*, "class")= chr [1:2] "bm_scale" "bru_mapper"
#>  - attr(*, "class")= chr [1:2] "bm_list" "list"
str(m[2])
#> List of 1
#>  $ B: list()
#>   ..- attr(*, "class")= chr [1:2] "bm_scale" "bru_mapper"
#>  - attr(*, "class")= chr [1:2] "bm_list" "list"
```
