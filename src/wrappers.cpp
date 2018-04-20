#include <Rcpp.h>
#include "ktau.hpp"
#include "hoeffd.hpp"
#include "prho.hpp"
#include "srho.hpp"
#include "indep_test.hpp"

//' @export
// [[Rcpp::export]]
double ktau_cpp(const std::vector<double>& x,
                const std::vector<double>& y,
                const std::vector<double>& weights)
{
    return ktau::ktau(x, y, weights);
}


//' @export
// [[Rcpp::export]]
double hoeffd_cpp(const std::vector<double>& x,
                  const std::vector<double>& y,
                  const std::vector<double>& weights)
{
    return hoeffd::hoeffd(x, y, weights);
}

//' @export
// [[Rcpp::export]]
double prho_cpp(const std::vector<double>& x,
                const std::vector<double>& y,
                const std::vector<double>& weights)
{
    return prho::prho(x, y, weights);
}

//' @export
// [[Rcpp::export]]
double srho_cpp(const std::vector<double>& x,
                const std::vector<double>& y,
                const std::vector<double>& weights)
{
    return srho::srho(x, y, weights);
}

//' @export
// [[Rcpp::export]]
std::vector<double> rank_scores_cpp(const std::vector<double>& x,
                                    const std::vector<double>& weights)
{
    return utils::rank_scores(x, weights);
}

//' @export
// [[Rcpp::export]]
std::vector<double> bivariate_rank_cpp(const std::vector<double>& x,
                                       const std::vector<double>& y,
                                       const std::vector<double>& weights)
{
    return hoeffd::bivariate_rank(x, y, weights);
}

//' @export
// [[Rcpp::export]]
double indep_test_asymptotic_cpp(const std::vector<double>& x,
                                 const std::vector<double>& y,
                                 std::string method,
                                 const std::vector<double>& weights)
{
    return indep_test::indep_test_asymptotic(x, y, method, weights);
}


