#include <Rcpp.h>
#include "wdm.hpp"

// [[Rcpp::export]]
double wdm_cpp(const std::vector<double>& x,
               const std::vector<double>& y,
               std::string method,
               const std::vector<double>& weights,
               bool remove_missing)
{
    return wdm::wdm(x, y, method, weights, remove_missing);
}

std::vector<double> convert_vec(const Rcpp::NumericVector& x)
{
    return Rcpp::as<std::vector<double>>(x);
}

// [[Rcpp::export]]
Rcpp::NumericMatrix wdm_mat_cpp(const Rcpp::NumericMatrix& x,
                                std::string method,
                                const std::vector<double>& weights)
{
    using namespace Rcpp;
    size_t d = x.ncol();
    if (d == 1)
        throw std::runtime_error("x must have at least 2 columns.");

    NumericMatrix ms(d, d);
    for (size_t i = 0; i < x.cols(); i++) {
        for (size_t j = i; j < x.cols(); j++) {
            if (j == i) {
                ms(j, i) = 1.0;
                continue;
            }
            ms(i, j) = wdm::wdm(convert_vec(x(_, i)),
                                convert_vec(x(_, j)),
                                method,
                                weights);
            ms(j, i) = ms(i, j);
        }
    }

    return ms;
}

// [[Rcpp::export]]
Rcpp::List indep_test_cpp(const std::vector<double>& x,
                          const std::vector<double>& y,
                          std::string method,
                          const std::vector<double>& weights,
                          bool remove_missing,
                          std::string alternative)
{
    wdm::Indep_test test(x, y, method, weights, remove_missing, alternative);
    return Rcpp::List::create(
        Rcpp::Named("estimate") = test.estimate(),
        Rcpp::Named("statistic") = test.statistic(),
        Rcpp::Named("p_value") = test.p_value(),
        Rcpp::Named("n_eff") = test.n_eff(),
        Rcpp::Named("method") = method,
        Rcpp::Named("alternative") = alternative
    );
}

// [[Rcpp::export]]
void test()
{
    // input vectors
    std::vector<double> x{1, 3, 2, 5, 3, 2, 20, 15};
    std::vector<double> y{2, 12, 4, 7, 8, 14, 17, 6};

    // weights
    std::vector<double> w{1, 1, 2, 2, 1, 0, 0.5, 0.3};

    std::cout <<
        "unweighted Kendall's tau: " << wdm::wdm(x, y, "kendall") << std::endl;
    std::cout <<
        "weighted Kendall's tau: " <<  wdm::wdm(x, y, "kendall", w) << std::endl;

    // weighted independence test
    wdm::Indep_test test(x, y, "kendall", w);
    std::cout << "statistic: " << test.statistic() << std::endl;
    std::cout << "p-value: " << test.p_value() << std::endl;
}
