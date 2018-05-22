
wdm
===

[![Travis build status](https://travis-ci.org/tnagler/wdm-r.svg?branch=master)](https://travis-ci.org/tnagler/wdm-r) [![AppVeyor build status](https://ci.appveyor.com/api/projects/status/github/tnagler/wdm-r?branch=master&svg=true)](https://ci.appveyor.com/project/tnagler/wdm-r)

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
#> [1] 0.08563521
wdm(x, y, method = "kendall", weights = w)  # weighted
#> [1] 0.1048728
```

#### dependence in a matrix

``` r
x <- matrix(rnorm(100 * 3), 100, 3)
wdm(x, method = "spearman")               # unweighted
#>            [,1]       [,2]       [,3]
#> [1,] 1.00000000 0.04116412 0.04162016
#> [2,] 0.04116412 1.00000000 0.06217822
#> [3,] 0.04162016 0.06217822 1.00000000
wdm(x, method = "spearman", weights = w)  # weighted
#>             [,1]        [,2]        [,3]
#> [1,]  1.00000000  0.02223874 -0.05844855
#> [2,]  0.02223874  1.00000000 -0.04757413
#> [3,] -0.05844855 -0.04757413  1.00000000
```

#### independence tests

``` r
x <- rnorm(100)
y <- rpois(100, 1)  # all but Hoeffding's D can handle ties
w <- runif(100)
indep_test(x, y, method = "kendall")               # unweighted
#>    estimate statistic   p_value n_eff  method alternative
#> 1 -0.102849 -1.199541 0.2303178   100 kendall   two-sided
indep_test(x, y, method = "kendall", weights = w)  # weighted
#>     estimate statistic    p_value    n_eff  method alternative
#> 1 -0.1700221 -1.701925 0.08876937 76.51569 kendall   two-sided
```
