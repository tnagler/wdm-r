
wdm
===

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
#> [1] -0.07064812
wdm(x, y, method = "kendall", weights = w)  # weighted
#> [1] -0.02955453
```

#### dependence in a matrix

``` r
x <- matrix(rnorm(100 * 3), 100, 3)
wdm(x, method = "spearman")               # unweighted
#>              [,1]       [,2]         [,3]
#> [1,]  1.000000000 0.24799280 -0.006552655
#> [2,]  0.247992799 1.00000000  0.036759676
#> [3,] -0.006552655 0.03675968  1.000000000
wdm(x, method = "spearman", weights = w)  # weighted
#>             [,1]      [,2]        [,3]
#> [1,]  1.00000000 0.3320806 -0.09031678
#> [2,]  0.33208061 1.0000000  0.03356930
#> [3,] -0.09031678 0.0335693  1.00000000
```

#### independence tests

``` r
x <- rnorm(100)
y <- rpois(100, 1)  # all but Hoeffding's D can handle ties
w <- runif(100)
indep_test(x, y, method = "kendall")               # unweighted
#>    estimate statistic   p_value n_eff  method alternative
#> 1 0.1166794  1.419237 0.1558298   100 kendall   two-sided
indep_test(x, y, method = "kendall", weights = w)  # weighted
#>    estimate statistic   p_value    n_eff  method alternative
#> 1 0.1089319  1.158999 0.2464567 73.47093 kendall   two-sided
```
