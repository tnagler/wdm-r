---
output: github_document
---


# wdm

[![R-CMD-check](https://github.com/tnagler/wdm-r/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/tnagler/wdm-r/actions/workflows/R-CMD-check.yaml)
[![CRAN status](https://www.r-pkg.org/badges/version/wdm)](https://cran.r-project.org/package=wdm)
![CRAN downloads](http://cranlogs.r-pkg.org/badges/wdm)

R interface to the [wdm](https://github.com/tnagler/wdm) C++ library, which
provides efficient implementations of weighted dependence measures and related
independence tests:

   * Pearson's rho
   * Spearman's rho
   * Kendall's tau
   * Blomqvist's beta
   * Hoeffding's D
   * Chatterjee's xi

All measures are computed in *O(n* log *n)* time, where *n* is the number of
observations.

For a detailed description of the functionality, see the
[API documentation](https://tnagler.github.io/wdm-r/).


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

This repo contains [wdm](https://github.com/tnagler/wdm) as a submodule. For
a full clone use
``` shell
git clone --recurse-submodules <repo-address>
```

### Examples


``` r
library(wdm)
```

#####  Dependence between two vectors

``` r
x <- rnorm(100)
y <- rpois(100, 1)
w <- runif(100)
wdm(x, y, method = "kendall")               # unweighted
#> [1] 0.04547414
wdm(x, y, method = "kendall", weights = w)  # weighted
#> [1] 0.07764891
```

#####  Dependence in a matrix

``` r
x <- matrix(rnorm(100 * 3), 100, 3)
wdm(x, method = "spearman")               # unweighted
#>             [,1]        [,2]        [,3]
#> [1,]  1.00000000  0.10045005 -0.03279928
#> [2,]  0.10045005  1.00000000 -0.02744674
#> [3,] -0.03279928 -0.02744674  1.00000000
wdm(x, method = "spearman", weights = w)  # weighted
#>             [,1]       [,2]        [,3]
#> [1,]  1.00000000 0.22080359 -0.05369922
#> [2,]  0.22080359 1.00000000  0.01192067
#> [3,] -0.05369922 0.01192067  1.00000000
```

Chatterjee's xi measures the dependence of the second argument on the first,
so reversing the arguments can change the result:


``` r
wdm(x[, 1], x[, 2], method = "chatterjee")
#> [1] 0.01740174
wdm(x[, 2], x[, 1], method = "chatterjee")
#> [1] 0.03510351
```

##### Independence test

``` r
x <- rnorm(100)
y <- rpois(100, 1)
w <- runif(100)
indep_test(x, y, method = "kendall")               # unweighted
#>     estimate statistic   p_value n_eff  method alternative
#> 1 0.02140979 0.2799245 0.7795354   100 kendall   two-sided
indep_test(x, y, method = "kendall", weights = w)  # weighted
#>     estimate statistic   p_value   n_eff  method alternative
#> 1 0.04616801 0.5015067 0.6160145 71.2503 kendall   two-sided
```

For Chatterjee's xi, `x` is the predictor and `y` is the response. Its natural
one-sided alternative is `"greater"`:


``` r
indep_test(x, y, method = "chatterjee", alternative = "greater")
#>      estimate statistic   p_value n_eff     method alternative
#> 1 -0.04496517 -0.576012 0.7176965   100 chatterjee     greater
```
