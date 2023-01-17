
# wdm

[![R-CMD-check](https://github.com/tnagler/wdm-r/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/tnagler/wdm-r/actions/workflows/R-CMD-check.yaml)
[![Coverage
status](https://codecov.io/gh/tnagler/wdm-r/branch/master/graph/badge.svg)](https://codecov.io/github/tnagler/wdm-r?branch=master)
[![CRAN
status](https://www.r-pkg.org/badges/version/wdm)](https://cran.r-project.org/package=wdm)
![CRAN downloads](http://cranlogs.r-pkg.org/badges/wdm)

R interface to the [wdm](https://github.com/tnagler/wdm) C++ library,
which provides efficient implementations of weighted dependence measures
and related independence tests:

- Pearsons’s rho
- Spearmans’s rho
- Kendall’s tau
- Blomqvist’s beta
- Hoeffding’s D

All measures are computed in *O(n* log *n)* time, where *n* is the
number of observations.

For a detailed description of the functionality, see the [API
documentation](https://tnagler.github.io/wdm-r/).

### Installation

- the stable release from CRAN:

``` r
install.packages("wdm")
```

- the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
install_submodule_git <- function(x, ...) {
  install_dir <- tempfile()
  system(paste("git clone --recursive", shQuote(x), shQuote(install_dir)))
  devtools::install(install_dir, ...)
}
install_submodule_git("https://github.com/tnagler/wdm-r")
```

### Cloning

This repo contains [wdm](https://github.com/tnagler/wdm) as a submodule.
For a full clone use

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
#> [1] -0.1206219
wdm(x, y, method = "kendall", weights = w)  # weighted
#> [1] -0.01751542
```

##### Dependence in a matrix

``` r
x <- matrix(rnorm(100 * 3), 100, 3)
wdm(x, method = "spearman")               # unweighted
#>            [,1]       [,2]       [,3]
#> [1,]  1.0000000 -0.1842904 -0.0060006
#> [2,] -0.1842904  1.0000000  0.1281728
#> [3,] -0.0060006  0.1281728  1.0000000
wdm(x, method = "spearman", weights = w)  # weighted
#>             [,1]       [,2]       [,3]
#> [1,]  1.00000000 -0.1435049 0.03738866
#> [2,] -0.14350490  1.0000000 0.23603074
#> [3,]  0.03738866  0.2360307 1.00000000
```

##### Independence test

``` r
x <- rnorm(100)
y <- rpois(100, 1)  # all but Hoeffding's D can handle ties
w <- runif(100)
indep_test(x, y, method = "kendall")               # unweighted
#>     estimate statistic   p_value n_eff  method alternative
#> 1 0.05391935 0.6092239 0.5423761   100 kendall   two-sided
indep_test(x, y, method = "kendall", weights = w)  # weighted
#>     estimate statistic   p_value    n_eff  method alternative
#> 1 0.01193407  0.115073 0.9083872 75.82538 kendall   two-sided
```
