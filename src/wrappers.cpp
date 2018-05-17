#include <Rcpp.h>
#include "wdm/wdm.hpp"

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


// [[Rcpp::export]]
double indeptest_cpp(const std::vector<double>& x,
                     const std::vector<double>& y,
                     std::string method,
                     const std::vector<double>& weights)
{
    return indeptest(x, y, method, weights);
}
