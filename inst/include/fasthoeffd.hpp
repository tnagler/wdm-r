#pragma once

#include "utils.hpp"

namespace fasthoeffd {


inline double fast_hoeffd(std::vector<double> x,
                          std::vector<double> y,
                          std::vector<double> weights = std::vector<double>())
{
    // check input sizes
    size_t n = x.size();
    if (y.size() != n)
        throw std::runtime_error("lengths of x and y must match.");
    bool weighted = (weights.size() > 0);
    if (weighted && (weights.size() != n))
        throw std::runtime_error("lengths of x, y, and weights must match.");

    if (weights.size() == 0)
        weights = std::vector<double>(n, 1);

    // 1. Compute the ranks
    std::vector<double> R_X = utils::compute_ranks(x, weights);
    std::vector<double> S_X = utils::compute_ranks(x, utils::pow(weights, 2));
    std::vector<double> R_Y = utils::compute_ranks(y, weights);
    std::vector<double> S_Y = utils::compute_ranks(y, utils::pow(weights, 2));

    // 2. Sort all vectors according to first variable, breaking ties with second.
    std::vector<double> R_tmp, w_tmp;
    utils::sort_all(R_X, R_Y, weights);
    w_tmp = std::vector<double>();
    utils::sort_all(S_X, S_Y, w_tmp);

    // 3. Compute (weighted) bivariate ranks (number of points w/ both columns
    // less than the ith row). Computetd by sorting the second variable again in
    // descending order while counting inversions.
    std::vector<double> R_XY(n, 0.0), S_XY(n, 0.0), Q_XY(n, 0.0), U_XY(n, 0.0);

    R_tmp = R_Y;
    w_tmp = weights;
    utils::merge_sort_count_per_element(R_tmp, R_X, w_tmp);
    R_tmp = R_Y;
    w_tmp = weights;
    utils::merge_sort_count_per_element(R_tmp, S_X, w_tmp);

    R_tmp = R_Y;
    w_tmp = utils::pow(weights, 4);
    utils::merge_sort_count_per_element(R_tmp, w_tmp, U_XY);

    R_tmp = R_Y;
    w_tmp = utils::pow(weights, 3);
    utils::merge_sort_count_per_element(R_tmp, w_tmp, Q_XY);

    w_tmp = utils::pow(weights, 2);
    utils::merge_sort_count_per_element(S_Y, w_tmp, S_XY);

    utils::merge_sort_count_per_element(R_Y, weights, R_XY);

    // 4. Compute (weighted) Hoeffdings' D
    double A = 0.0, B = 0.0, C = 0.0;
    for (size_t i = 0; i < n; i++) {
        A += (R_XY[i] * R_XY[i] - S_XY[i]) * weights[i];
        B += (
            (R_X[i] * R_Y[i] - S_XY[i]) * R_XY[i] -
            S_XY[i] * (R_X[i] + R_Y[i]) + 2 * Q_XY[i]
        ) * weights[i];
        C += (
            (R_X[i] * R_X[i] - S_X[i]) * (R_Y[i] * R_Y[i] - S_Y[i])  -
            4 * ((R_X[i] * R_Y[i] - S_XY[i]) * S_XY[i] -
            Q_XY[i] * (R_X[i] + R_Y[i]) + 2 * U_XY[i]) -
            2 * (S_XY[i] * S_XY[i] - U_XY[i])
        ) * weights[i];
    }

    double D = 0;
    D += A / (utils::perm_sum(weights, 3) * 6);
    D -= 2 * B / (utils::perm_sum(weights, 4) * 24);
    D += C / (utils::perm_sum(weights, 5) * 120);

    return 30.0 * D;
}

}
