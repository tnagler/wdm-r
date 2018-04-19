#include <RcppEigen.h>
#include <Eigen/Dense>
#include "fastktau.hpp"
#include "fasthoeffd.hpp"

//' @export
// [[Rcpp::export]]
double fast_ktau_cpp(const Eigen::MatrixXd& x, const Eigen::VectorXd& weights)
{
    size_t n = x.rows();
    std::vector<double> xx(n), yy(n), w(weights.size());
    Eigen::VectorXd::Map(&xx[0], n) = x.col(0);
    Eigen::VectorXd::Map(&yy[0], n) = x.col(1);
    Eigen::VectorXd::Map(&w[0], weights.size()) = weights;
    return fastktau::fast_ktau(xx, yy, w);
}


//' @export
// [[Rcpp::export]]
double fast_hoeffd_cpp(const Eigen::MatrixXd& x, const Eigen::VectorXd& weights)
{
    size_t n = x.rows();
    std::vector<double> xx(n), yy(n), w(weights.size());
    Eigen::VectorXd::Map(&xx[0], n) = x.col(0);
    Eigen::VectorXd::Map(&yy[0], n) = x.col(1);
    Eigen::VectorXd::Map(&w[0], weights.size()) = weights;
    return fasthoeffd::fast_hoeffd(xx, yy, w);
}


//' @export
// [[Rcpp::export]]
double fast_hoeffd2_cpp(const Eigen::MatrixXd& x, const Eigen::VectorXd& weights)
{
    //! calculates the pair-wise Hoeffding's D.
    size_t n = x.rows();

    // Compute the ranks
    Eigen::MatrixXd R = x;
    std::vector<size_t> order(n);
    std::vector<double> w(weights.size());

    for (size_t i = 0; i < n; i++)
        order[i] = i;
    std::sort(order.begin(), order.end(), [&](size_t i, size_t j) {
        return (x(i, 1) < x(j, 1));
    });
    for (size_t i = 0; i < n; i++)
        R(order[i], 1) = static_cast<double>(i + 1);

    for (size_t i = 0; i < n; i++)
        order[i] = i;
    std::sort(order.begin(), order.end(), [&](size_t i, size_t j) {
        return (x(i, 0) < x(j, 0));
    });
    for (size_t i = 0; i < n; i++)
        R(order[i], 0) = static_cast<double>(i + 1);

    // Compute Q, with Qi the number of points with both columns less than
    // their ith value
    Eigen::VectorXd Q(n);
    Eigen::MatrixXd tmp(n, 2);
    for (size_t i = 0; i < n; i++) {
        tmp.col(0) = Eigen::VectorXd::Constant(n, x(i, 0));
        tmp.col(1) = Eigen::VectorXd::Constant(n, x(i, 1));
        tmp = (x - tmp).unaryExpr([](double v) {
            double res = 0.0;
            if (v < 0.0) {
                res = 1.0;
            }
            return res;
        });
        Q(i) = tmp.rowwise().prod().sum();
    }

    Eigen::Matrix<double, Eigen::Dynamic, 2> ones = Eigen::MatrixXd::Ones(n, 2);
    double A = (R - ones).cwiseProduct(R - 2 * ones).rowwise().prod().sum();
    double B = (R - 2 * ones).rowwise().prod().cwiseProduct(Q).sum();
    double C = Q.cwiseProduct(Q - ones.col(0)).sum();

    double D = (A - 2 * (n - 2) * B + (n - 2) * (n - 3) * C);
    D /= (n * (n - 1) * (n - 2) * (n - 3) * (n - 4));

    return 30.0 * D;
}
