// Copyright © 2018 Thomas Nagler
//
// This file is part of the wdm library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory
// or https://github.com/tnagler/wdmcpp/blob/master/LICENSE.

#pragma once

#include "utils.hpp"

namespace wdm {

//! fast calculation of the weighted Kendall's tau.
//! @param x, y input data.
//! @param weights an optional vector of weights for the data.
inline double ktau(std::vector<double> x,
                   std::vector<double> y,
                   std::vector<double> weights = std::vector<double>())
{
    wdm_utils::check_sizes(x, y, weights);

    // 1.1 Sort x, y, and weights in x order; break ties in according to y.
    wdm_utils::sort_all(x, y, weights);

    // 1.2 Count pairs of tied x and simultaneous ties in x and y.
    double ties_x = wdm_utils::count_ties(x, weights);
    double ties_both = wdm_utils::count_joint_ties(x, y, weights);

    // 2.1 Sort y again and count exchanges (= number of discordant pairs).
    double num_d = 0.0;
    wdm_utils::merge_sort(y, weights, num_d);

    // 2.2 Count pairs of tied y.
    double ties_y = wdm_utils::count_ties(y, weights);

    // 3. Calculate Kendall's tau.
    if (weights.size() == 0)
        weights = std::vector<double>(x.size(), 1.0);
    double num_pairs = wdm_utils::perm_sum(weights, 2);
    double num_c = num_pairs - (num_d + ties_x + ties_y - ties_both);
    double tau = num_c - num_d;
    tau /= std::sqrt((num_pairs - ties_x) * (num_pairs - ties_y));

    return tau;
}

}
