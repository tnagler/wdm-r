
wdm
===

[![Travis build status](https://travis-ci.org/tnagler/wdm-r.svg?branch=master)](https://travis-ci.org/tnagler/wdm-r) [![AppVeyor build status](https://ci.appveyor.com/api/projects/status/github/tnagler/wdm-r?branch=master&svg=true)](https://ci.appveyor.com/project/tnagler/wdm-r) [![Coverage status](https://codecov.io/gh/tnagler/wdm-r/branch/master/graph/badge.svg)](https://codecov.io/github/tnagler/wdm-r?branch=master) [![CRAN status](https://www.r-pkg.org/badges/version/wdm)](https://cran.r-project.org/package=wdm)

R interface to the [wdm](https://github.com/tnagler/wdm) C++ library, which provides efficient implementations of weighted dependence measures and related independence tests:

-   Pearsons's rho
-   Spearmans's rho
-   Kendall's tau
-   Blomqvist's beta
-   Hoeffding's D

All measures are computed in *O(n* log *n)* time, where *n* is the number of observations. For a detailed description of the functionality, see the [API documentation](https://tnagler.github.io/wdm-r/).

### Installation

You can install the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
install_submodule_git <- function(x, ...) {
  install_dir <- tempfile()
  system(paste("git clone --recursive", shQuote(x), shQuote(install_dir)))
  devtools::install(install_dir, ...)
}
install_submodule_git("https://github.com/tnagler/wdm-r")
```

This assumes you have [git](https://git-scm.com) client installed locally.

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
#> [1] 0.03877453
wdm(x, y, method = "kendall", weights = w)  # weighted
#> [1] 0.02834627
```

##### Dependence in a matrix

``` r
x <- matrix(rnorm(100 * 3), 100, 3)
wdm(x, method = "spearman")               # unweighted
#>             [,1]       [,2]        [,3]
#> [1,]  1.00000000 0.09048905 -0.10154215
#> [2,]  0.09048905 1.00000000  0.04774077
#> [3,] -0.10154215 0.04774077  1.00000000
wdm(x, method = "spearman", weights = w)  # weighted
#>             [,1]       [,2]        [,3]
#> [1,]  1.00000000 0.07679064 -0.01569089
#> [2,]  0.07679064 1.00000000  0.03786568
#> [3,] -0.01569089 0.03786568  1.00000000
```

##### Independence test

``` r
x <- rnorm(100)
y <- rpois(100, 1)  # all but Hoeffding's D can handle ties
w <- runif(100)
indep_test(x, y, method = "kendall")               # unweighted
#>     estimate statistic    p_value n_eff  method alternative
#> 1 -0.1477247 -1.756304 0.07903658   100 kendall   two-sided
indep_test(x, y, method = "kendall", weights = w)  # weighted
#>     estimate statistic   p_value    n_eff  method alternative
#> 1 -0.1352065 -1.357386 0.1746587 73.98258 kendall   two-sided
```
