# wdm 0.3.0

* update the bundled C++ backend to wdm 0.3.0.
* add weighted Chatterjee's xi and its asymptotic independence test.
* make predictor-tie handling for Chatterjee's xi reproducible by default and
  expose explicit seeds for alternative tie orderings.
* document that Hoeffding's D does not support tied observations.

# wdm 0.2.4

* add function `rank_wtd()` for computing weighted ranks.
 
 # wdm 0.2.3

* avoid bitwise operations on Boolean variables. 
 
 
# wdm 0.2.2

* fix computation of average weighted rank (no longer dependent on input order).
  Primarily affects method "spearman".


# wdm 0.2.1

* keep debugging symbols after package install (following a request by Prof. Ripley).


# wdm 0.2.0

* use more intuitive formula for Pearson's rho that is linear in weights.


# wdm 0.1.2

* fix typos and formatting issues in manual pages.


# wdm 0.1.1

BUG FIX

* inline `normalCDF()` to avoid multiple definition errors when `LinkingTo` wdm.


# wdm 0.1.0

* Initial release.
