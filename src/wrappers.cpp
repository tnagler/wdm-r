#include <Rcpp.h>
#include "wdmcpp/wdm.hpp"

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
double wdm_cpp(const std::vector<double>& x,
               const std::vector<double>& y,
               std::string method,
               const std::vector<double>& weights)
{
    return wdm::wdm(x, y, method, weights);
}


//' @export
// [[Rcpp::export]]
double indeptest_cpp(const std::vector<double>& x,
                     const std::vector<double>& y,
                     std::string method,
                     const std::vector<double>& weights)
{
    return indeptest(x, y, method, weights);
}

//' @export
// [[Rcpp::export]]
void test()
{
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
}
