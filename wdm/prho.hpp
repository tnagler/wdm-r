// Copyright © 2018 Thomas Nagler
//
// This file is part of the wdm library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory
// or https://github.com/tnagler/wdmcpp/blob/master/LICENSE.

#pragma once

#include "utils.hpp"


namespace wdm {

//! fast calculation of the weighted Pearson's correlation.
//! @param x, y input data.
//! @param weights an optional vector of weights for the data.
inline double prho(std::vector<double> x,
                   std::vector<double> y,
                   std::vector<double> weights = std::vector<double>())
{
    wdm_utils::check_sizes(x, y, weights);
    size_t n = x.size();
    if (weights.size() == 0)
        weights = std::vector<double>(x.size(), 1.0);

    // calculate means of x and y
    double mu_x = 0.0, mu_y = 0.0, w_sum = 0.0;
    for (size_t i = 0; i < n; i++) {
        mu_x += x[i] * weights[i];
        mu_y += y[i] * weights[i];
        w_sum += weights[i];
    }
    mu_x /= w_sum;
    mu_y /= w_sum;

    // substract mean from x and y
    for (size_t i = 0; i < n; i++) {
        x[i] -= mu_x;
        y[i] -= mu_y;
    }

    // compute variances and covariance
    double v_x = 0.0, v_y = 0.0, cov = 0.0, wi_sq;
    for (size_t i = 0; i < n; i++) {
        wi_sq = weights[i] * weights[i];
        v_x += x[i] * x[i] * wi_sq;
        v_y += y[i] * y[i] * wi_sq;
        cov += x[i] * y[i] * wi_sq;
    }

    // compute correlation
    return cov / std::sqrt(v_x * v_y);
}

}
