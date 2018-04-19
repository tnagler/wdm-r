#pragma once

#include "utils.hpp"

namespace fastktau {

//! fast calculation of the weighted Kendall's tau.
//! @param x, y input data.
//! @param weights an optional vector of weights for the data.
inline double fast_ktau(
        std::vector<double> x,
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

    // 0.1 Normalize weights such that they sum to choose(n, 2).
    utils::normalize_weights(weights);

    // 1.1 Sort x, y, and weights in x order; break ties in according to y.
    utils::sort_all(x, y, weights);

    // 1.2 Count pairs of tied x and simultaneous ties in x and y.
    double ties_x = utils::count_ties(x, weights);
    double ties_both = utils::count_joint_ties(x, y, weights);

    // 2.1 Sort y again and count exchanges (= number of discordant pairs).
    double num_d = 0.0;
    utils::merge_sort(y, weights, num_d);

    // 2.2 Count pairs of tied y.
    double ties_y = utils::count_ties(y, weights);

    // 3. Calculate Kendall's tau.
    double num_pairs = 0.5 * n * (n - 1);
    double num_c = num_pairs - (num_d + ties_x + ties_y - ties_both);
    double tau = num_c - num_d;
    tau /= std::sqrt((num_pairs - ties_x) * (num_pairs - ties_y));

    return tau;
}

}

