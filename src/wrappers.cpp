#include <Rcpp.h>
#include "ktau.hpp"
#include "hoeffd.hpp"

//' @export
// [[Rcpp::export]]
double fast_ktau_cpp(const std::vector<double>& x,
                     const std::vector<double>& y,
                     const std::vector<double>& weights)
{
    return ktau::ktau(x, y, weights);
}


//' @export
// [[Rcpp::export]]
double fast_hoeffd_cpp(const std::vector<double>& x,
                       const std::vector<double>& y,
                       const std::vector<double>& weights)
{
    return hoeffd::hoeffd(x, y, weights);
}
