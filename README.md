# wdm
A header-only C++ library for weighted dependence measures

Provides an efficient implementation of weighted dependence measures:

   * Pearsons's rho: `prho()`
   * Spearmans's rho: `srho()`
   * Kendall's tau: `ktau()`
   * Blomqvist's beta: `bbeta()`
   * Hoeffding's D: `hoeffd()`

All functions take at least two arguments for the two input vectors to compute 
the dependence measure on, plus an optional argument to supply weights for 
each observation. Additionally asymptotic p-values for all measures can be 
calculated with `indep_test()`.

### Example

``` cpp
#include "wdm.hpp"

// input vectors
std::vector<double> x{1, 3, 2, 5, 3, 2, 20, 15};
std::vector<double> y{2, 12, 4, 7, 8, 14, 17, 6};

// weights
std::vector<double> w{1, 1, 2, 2, 1, 0, 0.5, 0.3};

// unweighted Kendall's tau
std::cout << "Kendall's tau: " << wdm::ktau(x, y) << std::endl;
// or: wdm::wdm(x, y, "ktau")
std::cout << "p-value: " << wdm::indeptest(x, y, "ktau") << std::endl;

// weighted Kendall's tau
std::cout << "Kendall's tau: " << wdm::ktau(x, y, w) << std::endl;
// or: wdm::wdm(x, y, "ktau", w)
std::cout << "p-value: " << wdm::indeptest(x, y, "ktau", w) << std::endl;
```
Output:
```
Kendall's tau: 0.2965
p-value: 0.208413
Kendall's tau: 0.550633
p-value: 0.0557333
```
