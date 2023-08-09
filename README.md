
# wdm

[![R-CMD-check](https://github.com/tnagler/wdm-r/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/tnagler/wdm-r/actions/workflows/R-CMD-check.yaml)
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
#> [1] -0.03093257
wdm(x, y, method = "kendall", weights = w)  # weighted
#> [1] 0.04835766
```

##### Dependence in a matrix

``` r
x <- matrix(rnorm(100 * 3), 100, 3)
wdm(x, method = "spearman")               # unweighted
#>             [,1]      [,2]        [,3]
#> [1,]  1.00000000 0.2194659 -0.05435344
#> [2,]  0.21946595 1.0000000  0.11401140
#> [3,] -0.05435344 0.1140114  1.00000000
wdm(x, method = "spearman", weights = w)  # weighted
#>            [,1]      [,2]       [,3]
#> [1,]  1.0000000 0.2575236 -0.1689466
#> [2,]  0.2575236 1.0000000  0.1197442
#> [3,] -0.1689466 0.1197442  1.0000000
```

##### Independence test

``` r
x <- rnorm(100)
y <- rpois(100, 1)  # all but Hoeffding's D can handle ties
w <- runif(100)
indep_test(x, y, method = "kendall")               # unweighted
#>    estimate statistic   p_value n_eff  method alternative
#> 1 0.1278922  1.532215 0.1254693   100 kendall   two-sided
indep_test(x, y, method = "kendall", weights = w)  # weighted
#>    estimate statistic    p_value   n_eff  method alternative
#> 1 0.1704296  1.779486 0.07516007 79.6939 kendall   two-sided
```
