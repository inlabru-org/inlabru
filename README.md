
<!-- README.md is generated from README.Rmd. Please edit that file -->
inlabru
=======

[![Build Status](https://travis-ci.org/fbachl/inlabru.svg?branch=devel)](https://travis-ci.org/fbachl/inlabru) [![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/inlabru)](https://cran.r-project.org/package=inlabru)

The goal of [inlabru](http://inlabru.org) is to facilitate spatial modeling using integrated nested Laplace approximation via the [R-INLA package](http://www.r-inla.org). Additionally, implements a log Gaussian Cox process likelihood for modeling univariate and spatial point processes based on ecological survey data. See Yuan Yuan, Fabian E. Bachl, Finn Lindgren, David L. Borchers, Janine B. Illian, Stephen T. Buckland, Havard Rue, Tim Gerrodette (2017), [arXiv](https://arxiv.org/abs/1604.06013).

Installation
------------

You can install the current [CRAN](https://CRAN.R-project.org) version of inlabru:

``` r
install.packages("inlabru")
```

You can install the latest bugfix release of inlabru from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("fbachl/inlabru", ref="master")
```

You can install the development version of inlabru from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("fbachl/inlabru", ref="devel")
```

Example
-------

This is a basic example which shows you how to solve a common problem:

``` r
## basic example code
```
