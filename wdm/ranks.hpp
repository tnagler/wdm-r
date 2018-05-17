// Copyright © 2018 Thomas Nagler
//
// This file is part of the wdm library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory
// or https://github.com/tnagler/wdmcpp/blob/master/LICENSE.

#pragma once 

#include "utils.hpp"

namespace wdm {
//! computes ranks (such that smallest element has rank 0), assigning average
//! ranks for ties.
//! @param x input vector.
//! @param ties_method `"min"` (default) assigns all tied values the minimum
//!   score; `"average"` assigns the average score.
//! @param weights (optional), weights for each observation.
//! @return a vector containing the ranks of each element in `x`.
std::vector<double> rank_scores(
    std::vector<double> x,
    std::vector<double> weights = std::vector<double>(),
    std::string ties_method = "min")
{
    size_t n = x.size();
    if (weights.size() == 0)
        weights = std::vector<double>(n, 1.0);

    std::vector<size_t> perm = wdm_utils::get_order(x);

    double w_acc = 0.0, w_batch;
    if ((ties_method != "min") && (ties_method != "average"))
        throw std::runtime_error("ties_method must be either 'min' or 'average.");
    for (size_t i = 0, reps; i < n; i += reps) {
        // find replications
        reps = 1;
        w_batch = 0.0;
        while ((i + reps < n) && (x[perm[i]] == x[perm[i + reps]])) {
            w_batch += weights[perm[i]];
            ++reps;
        }

        // assign average rank of the tied values
        for (size_t k = 0; k < reps; ++k) {
            if (ties_method == "min")
                x[perm[i + k]] = w_acc;
            else
                x[perm[i + k]] = w_acc + w_batch / 2.0;
        }
        
        // accumulate weights for current batch
        for (size_t k = 0; k < reps; ++k)
            w_acc += weights[perm[i + k]];
    }

    return x;
}

//! computes the bivariate rank of a pair of vectors (starting at 0).
//! @param x first input vector.
//! @param y second input vecotr.
//! @param weights (optional), weights for each observation.
inline std::vector<double> bivariate_rank(
        std::vector<double> x,
        std::vector<double> y,
        std::vector<double> weights = std::vector<double>())
{
    wdm_utils::check_sizes(x, y, weights);

    // get inverse of permutation that brings x in ascending order
    std::vector<size_t> perm_x = wdm_utils::get_order(x);
    perm_x = wdm_utils::invert_permutation(perm_x);

    // sort x, y, and weights according to x, breaking ties with y
    wdm_utils::sort_all(x, y, weights);

    // get inverse of permutation that brings y in descending order
    std::vector<size_t> perm_y = wdm_utils::get_order(y, false);
    perm_y = wdm_utils::invert_permutation(perm_y);

    // sort y in descending order counting inversions
    std::vector<double> counts(y.size(), 0.0);
    wdm_utils::merge_sort_count_per_element(y, weights, counts);

    // bring counts back in original order
    std::vector<double> counts_tmp = counts;
    for (size_t i = 0; i < counts.size(); i++)
        counts[i] = counts_tmp[perm_y[perm_x[i]]];

    return counts;
}
    
}
