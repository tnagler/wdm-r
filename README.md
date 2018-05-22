
wdm
===

[![Travis build status](https://travis-ci.org/tnagler/wdm-r.svg?branch=master)](https://travis-ci.org/tnagler/wdm-r) [![AppVeyor build status](https://ci.appveyor.com/api/projects/status/github/tnagler/wdm-r?branch=master&svg=true)](https://ci.appveyor.com/project/tnagler/wdm-r) [![Coverage status](https://codecov.io/gh/tnagler/wdm-r/branch/master/graph/badge.svg)](https://codecov.io/github/tnagler/wdm-r?branch=master)

R interface to the [wdm](https://github.com/tnagler/wdm) C++ library, which provides efficient implementations of weighted dependence measures and related independence tests:

-   Pearsons's rho
-   Spearmans's rho
-   Kendall's tau
-   Blomqvist's beta
-   Hoeffding's D

All measures are computed in *O(n* log *n)* time, where *n* is the number of observations.

Installation
------------

You can install the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("tnagler/wdm-r")
```

Examples
--------

``` r
library(wdm)
```

#### dependence between two vectors

``` r
x <- rnorm(100)
y <- rpois(100, 1)  # all but Hoeffding's D can handle ties
w <- runif(100)
wdm(x, y, method = "kendall")               # unweighted
#> [1] -0.07682849
wdm(x, y, method = "kendall", weights = w)  # weighted
#> [1] -0.08959182
```

#### dependence in a matrix

``` r
x <- matrix(rnorm(100 * 3), 100, 3)
wdm(x, method = "spearman")               # unweighted
#>             [,1]        [,2]        [,3]
#> [1,]  1.00000000 -0.07431143 -0.03681968
#> [2,] -0.07431143  1.00000000 -0.08433243
#> [3,] -0.03681968 -0.08433243  1.00000000
wdm(x, method = "spearman", weights = w)  # weighted
#>             [,1]        [,2]        [,3]
#> [1,]  1.00000000 -0.26920202 -0.06345475
#> [2,] -0.26920202  1.00000000 -0.06665639
#> [3,] -0.06345475 -0.06665639  1.00000000
```

#### independence tests

``` r
x <- rnorm(100)
y <- rpois(100, 1)  # all but Hoeffding's D can handle ties
w <- runif(100)
indep_test(x, y, method = "kendall")               # unweighted
#>    estimate statistic   p_value n_eff  method alternative
#> 1 0.0366728 0.4311852 0.6663337   100 kendall   two-sided
indep_test(x, y, method = "kendall", weights = w)  # weighted
#>     estimate statistic   p_value    n_eff  method alternative
#> 1 0.02995828 0.3055508 0.7599467 74.16138 kendall   two-sided
```
