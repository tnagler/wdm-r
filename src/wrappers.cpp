#define BOOST_MATH_PROMOTE_DOUBLE_POLICY false
#include <Rcpp.h>
#include "wdm.hpp"

// [[Rcpp::export]]
double wdm_cpp(const std::vector<double>& x,
               const std::vector<double>& y,
               std::string method,
               const std::vector<double>& weights,
               bool remove_missing,
               const std::vector<int>& seeds)
{
    return wdm::wdm(x, y, method, weights, remove_missing, seeds);
}

std::vector<double> convert_vec(const Rcpp::NumericVector& x)
{
    return Rcpp::as<std::vector<double>>(x);
}

// [[Rcpp::export]]
Rcpp::NumericMatrix wdm_mat_cpp(const Rcpp::NumericMatrix& x,
                                std::string method,
                                const std::vector<double>& weights,
                                bool remove_missing,
                                const std::vector<int>& seeds)
{
    using namespace Rcpp;
    int d = x.ncol();

    NumericMatrix ms(d, d);
    for (int i = 0; i < x.cols(); i++) {
        for (int j = i; j < x.cols(); j++) {
            ms(i, j) = wdm::wdm(convert_vec(x(_, i)),
                                convert_vec(x(_, j)),
                                method,
                                weights,
                                remove_missing,
                                seeds);
            if (j == i)
                continue;
            if (method == "chatterjee") {
                ms(j, i) = wdm::wdm(convert_vec(x(_, j)),
                                    convert_vec(x(_, i)),
                                    method,
                                    weights,
                                    remove_missing,
                                    seeds);
            } else {
                ms(j, i) = ms(i, j);
            }
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
                          std::string alternative,
                          const std::vector<int>& seeds,
                          bool y_continuous)
{
    wdm::Indep_test test(x,
                         y,
                         method,
                         weights,
                         remove_missing,
                         alternative,
                         seeds,
                         y_continuous);
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
std::vector<double> rank_wtd_cpp(
        std::vector<double> x,
        std::vector<double> weights,
        std::string ties_method = "min")
{
    return wdm::impl::rank(x, weights, ties_method);
}

// [[Rcpp::export]]
double perm_sum_cpp(const std::vector<double>& x, size_t k) {
    return wdm::utils::perm_sum(x, k);
}
