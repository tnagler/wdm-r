#include <Rcpp.h>
#include "wdm.hpp"

using namespace wdm;

//' @export
// [[Rcpp::export]]
double ktau_cpp(const std::vector<double>& x,
                const std::vector<double>& y,
                const std::vector<double>& weights)
{
    return ktau(x, y, weights);
}


//' @export
// [[Rcpp::export]]
double hoeffd_cpp(const std::vector<double>& x,
                  const std::vector<double>& y,
                  const std::vector<double>& weights)
{
    return hoeffd(x, y, weights);
}

//' @export
// [[Rcpp::export]]
double prho_cpp(const std::vector<double>& x,
                const std::vector<double>& y,
                const std::vector<double>& weights)
{
    return prho(x, y, weights);
}

//' @export
// [[Rcpp::export]]
double srho_cpp(const std::vector<double>& x,
                const std::vector<double>& y,
                const std::vector<double>& weights)
{
    return srho(x, y, weights);
}


//' @export
// [[Rcpp::export]]
double bbeta_cpp(const std::vector<double>& x,
                 const std::vector<double>& y,
                 const std::vector<double>& weights)
{
    return bbeta(x, y, weights);
}


//' @export
// [[Rcpp::export]]
std::vector<double> rank_scores_cpp(const std::vector<double>& x,
                                    const std::vector<double>& weights)
{
    return rank_scores(x, weights);
}

//' @export
// [[Rcpp::export]]
std::vector<double> bivariate_rank_cpp(const std::vector<double>& x,
                                       const std::vector<double>& y,
                                       const std::vector<double>& weights)
{
    return bivariate_rank(x, y, weights);
}

//' @export
// [[Rcpp::export]]
double indep_test_cpp(const std::vector<double>& x,
                                 const std::vector<double>& y,
                                 std::string method,
                                 const std::vector<double>& weights)
{
    return indep_test(x, y, method, weights);
}
