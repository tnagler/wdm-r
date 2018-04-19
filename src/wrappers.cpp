#include <Rcpp.h>
#include "ktau.hpp"
#include "hoeffd.hpp"
#include "prho.hpp"

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
