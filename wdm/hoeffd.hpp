// Copyright © 2018 Thomas Nagler
//
// This file is part of the wdm library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory
// or https://github.com/tnagler/wdmcpp/blob/master/LICENSE.

#pragma once

#include "utils.hpp"
#include "ranks.hpp"

namespace wdm {

//! fast calculation of the weighted Hoeffdings's D.
//! @param x, y input data.
//! @param weights an optional vector of weights for the data.
inline double hoeffd(std::vector<double> x,
                     std::vector<double> y,
                     std::vector<double> weights = std::vector<double>())
{
    if (x.size() < 5)
        throw std::runtime_error("number of observations must be at least 5.");
    wdm_utils::check_sizes(x, y, weights);

    // 1. Compute (weighted) ranks
    std::vector<double> R_X = rank_scores(x, weights);
    std::vector<double> R_Y = rank_scores(y, weights);
    std::vector<double> S_X, S_Y;
    if (weights.size() > 0) {
        S_X = rank_scores(x, wdm_utils::pow(weights, 2));
        S_Y = rank_scores(y, wdm_utils::pow(weights, 2));
    } else {
        S_X = R_X;
        S_Y = R_Y;
    }

    // 2. Compute (weighted) bivariate ranks (number of points w/ both columns
    // less than the ith row).
    std::vector<double> R_XY, S_XY, T_XY, U_XY;
    R_XY = bivariate_rank(x, y, weights);
    if (weights.size() > 0) {
        S_XY = bivariate_rank(x, y, wdm_utils::pow(weights, 2));
        T_XY = bivariate_rank(x, y, wdm_utils::pow(weights, 3));
        U_XY = bivariate_rank(x, y, wdm_utils::pow(weights, 4));
    } else {
        S_XY = R_XY;
        T_XY = R_XY;
        U_XY = R_XY;
    }


    // 3. Compute (weighted) Hoeffdings' D
    if (weights.size() == 0)
        weights = std::vector<double>(x.size(), 1.0);
    double A_1 = 0.0, A_2 = 0.0, A_3 = 0.0;
    for (size_t i = 0; i < x.size(); i++) {
        A_1 += (R_XY[i] * R_XY[i] - S_XY[i]) * weights[i];
        A_2 += (
            (R_X[i] * R_Y[i] - S_XY[i]) * R_XY[i] -
                S_XY[i] * (R_X[i] + R_Y[i]) + 2 * T_XY[i]
        ) * weights[i];
        A_3 += (
            (R_X[i] * R_X[i] - S_X[i]) * (R_Y[i] * R_Y[i] - S_Y[i])  -
                4 * ((R_X[i] * R_Y[i] - S_XY[i]) * S_XY[i] -
                T_XY[i] * (R_X[i] + R_Y[i]) + 2 * U_XY[i]) -
                2 * (S_XY[i] * S_XY[i] - U_XY[i])
        ) * weights[i];
    }
    double D = 0.0;
    D += A_1 / (wdm_utils::perm_sum(weights, 3) * 6);
    D -= 2 * A_2 / (wdm_utils::perm_sum(weights, 4) * 24);
    D += A_3 / (wdm_utils::perm_sum(weights, 5) * 120);

    return 30.0 * D;
}

}
