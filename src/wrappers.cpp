#include <Rcpp.h>
#include "wdm.hpp"

using namespace wdm;

// [[Rcpp::export]]
double ktau_cpp(const std::vector<double>& x,
                const std::vector<double>& y,
                const std::vector<double>& weights)
{
    return ktau(x, y, weights);
}


// [[Rcpp::export]]
double hoeffd_cpp(const std::vector<double>& x,
                  const std::vector<double>& y,
                  const std::vector<double>& weights)
{
    return hoeffd(x, y, weights);
}

// [[Rcpp::export]]
double prho_cpp(const std::vector<double>& x,
                const std::vector<double>& y,
                const std::vector<double>& weights)
{
    return prho(x, y, weights);
}

// [[Rcpp::export]]
double srho_cpp(const std::vector<double>& x,
                const std::vector<double>& y,
                const std::vector<double>& weights)
{
    return srho(x, y, weights);
}


// [[Rcpp::export]]
double bbeta_cpp(const std::vector<double>& x,
                 const std::vector<double>& y,
                 const std::vector<double>& weights)
{
    return bbeta(x, y, weights);
}


// [[Rcpp::export]]
std::vector<double> rank_scores_cpp(const std::vector<double>& x,
                                    const std::vector<double>& weights)
{
    return rank_scores(x, weights);
}

// [[Rcpp::export]]
std::vector<double> bivariate_rank_cpp(const std::vector<double>& x,
                                       const std::vector<double>& y,
                                       const std::vector<double>& weights)
{
    return bivariate_rank(x, y, weights);
}

// [[Rcpp::export]]
double wdm_cpp(const std::vector<double>& x,
               const std::vector<double>& y,
               std::string method,
               const std::vector<double>& weights)
{
    return wdm::wdm(x, y, method, weights);
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
double indeptest_cpp(const std::vector<double>& x,
                     const std::vector<double>& y,
                     std::string method,
                     const std::vector<double>& weights)
{
    return indeptest(x, y, method, weights);
}

