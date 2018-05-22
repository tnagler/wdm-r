
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

### Installation

You can install the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("tnagler/wdm-r")
```

### Cloning

This repo contains [wdm](https://github.com/tnagler/wdm) as a submodule. For a full clone use

``` shell
git clone --recurse-submodules <repo-address>
```

### Examples

``` r
library(wdm)
```

##### Dependence between two vectors

``` r
x <- rnorm(100)
y <- rpois(100, 1)  # all but Hoeffding's D can handle ties
w <- runif(100)
wdm(x, y, method = "kendall")               # unweighted
#> [1] 0.01850718
wdm(x, y, method = "kendall", weights = w)  # weighted
#> [1] 0.03520163
```

##### Dependence in a matrix

``` r
x <- matrix(rnorm(100 * 3), 100, 3)
wdm(x, method = "spearman")               # unweighted
#>             [,1]        [,2]        [,3]
#> [1,]  1.00000000 -0.09494149 -0.03687969
#> [2,] -0.09494149  1.00000000  0.05993399
#> [3,] -0.03687969  0.05993399  1.00000000
wdm(x, method = "spearman", weights = w)  # weighted
#>            [,1]        [,2]        [,3]
#> [1,]  1.0000000 -0.13994018 -0.13932619
#> [2,] -0.1399402  1.00000000  0.04920561
#> [3,] -0.1393262  0.04920561  1.00000000
```

##### Independence test

``` r
x <- rnorm(100)
y <- rpois(100, 1)  # all but Hoeffding's D can handle ties
w <- runif(100)
indep_test(x, y, method = "kendall")               # unweighted
#>      estimate  statistic   p_value n_eff  method alternative
#> 1 -0.01464269 -0.1711457 0.8641092   100 kendall   two-sided
indep_test(x, y, method = "kendall", weights = w)  # weighted
#>      estimate  statistic   p_value    n_eff  method alternative
#> 1 -0.01042896 -0.1039409 0.9172163 74.20691 kendall   two-sided
```
