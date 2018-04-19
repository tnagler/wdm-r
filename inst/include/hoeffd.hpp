#pragma once

#include "utils.hpp"

namespace hoeffd {

//! merge operation for a pair of vectors, counting inversions per element.
//! @param vec container for the sorted elements.
//! @param vec1, vec2 sorted input vectors to be merged.
//! @param weights container for the weights corresponding to sorted elements
//!   in `vec`; can be empty for unweighted counts.
//! @param weights1, weights2 weights corresponding to input vectors`vec1`,
//!   `vec2`; can be empty for unweighted counts.
//! @param counts container for the counts corresponding to sorted elements
//!   in `vec`.
//! @param counts1, counts2 counts corresponding to input vectors`vec1`,
//!   `vec2` to which (weighted) counts are added.
inline void merge_count_per_element(std::vector<double>& vec,
                                    const std::vector<double>& vec1,
                                    const std::vector<double>& vec2,
                                    std::vector<double>& weights,
                                    const std::vector<double>& weights1,
                                    const std::vector<double>& weights2,
                                    std::vector<double>& counts,
                                    const std::vector<double>& counts1,
                                    const std::vector<double>& counts2)
{
    double w_acc = 0.0;
    bool weighted = (weights.size() > 0);
    double w1_sum = 0.0;
    if (weighted) {
        for (size_t i = 0; i < weights1.size(); i++)
            w1_sum += weights1[i];
    }
    size_t i, j, k;
    for (i = 0, j = 0, k = 0; i < vec1.size() && j < vec2.size(); k++) {
        if (vec1[i] > vec2[j]) {
            vec[k] = vec1[i];
            counts[k] = counts1[i];
            if (weighted) {
                weights[k] = weights1[i];
                w_acc += weights1[i];
            }
            i++;
        } else {
            vec[k] = vec2[j];
            if (weighted) {
                counts[k] = counts2[j] + w1_sum - w_acc;
                weights[k] = weights2[j];
            } else {
                counts[k] = counts2[j] + vec1.size() - i;
            }
            j++;
        }
    }

    while (i < vec1.size()) {
        vec[k] = vec1[i];
        if (weighted)
            weights[k] = weights1[i];
        counts[k] = counts1[i];
        k++;
        i++;
    }

    while (j < vec2.size()) {
        vec[k] = vec2[j];
        if (weighted)
            weights[k] = weights2[j];
        counts[k] = counts2[j];
        k++;
        j++;
    }
}


//! sorting elements in a vector while counting inversions per element.
//! @param vec the vector to be sorted.
//! @param counts vector of counters to which the (weighted) number of inversions
//!   (per element) are added.
//! @param weights vector of weights corresponding to `vec`; can be empty for
//!   unweighted counts.
inline void merge_sort_count_per_element(std::vector<double>& vec,
                                         std::vector<double>& weights,
                                         std::vector<double>& counts)
{
    if (vec.size() > 1) {
        size_t n = vec.size();
        std::vector<double> vec1(vec.begin(), vec.begin() + n / 2);
        std::vector<double> vec2(vec.begin() + n / 2, vec.end());

        n = weights.size();
        std::vector<double> weights1(weights.begin(), weights.begin() + n / 2);
        std::vector<double> weights2(weights.begin() + n / 2, weights.end());

        n = counts.size();
        std::vector<double> counts1(counts.begin(), counts.begin() + n / 2);
        std::vector<double> counts2(counts.begin() + n / 2, counts.end());

        merge_sort_count_per_element(vec1, weights1, counts1);
        merge_sort_count_per_element(vec2, weights2, counts2);
        merge_count_per_element(vec, vec1, vec2,
                                weights, weights1, weights2,
                                counts, counts1, counts2);
    }
}

//! computes the bivariate rank of a pair of vectors (starting at 0).
//! @param x first input vector.
//! @param y second input vecotr.
//! @param (optional), weights for each observation.
//! @param return
inline std::vector<double> bivariate_rank(
        std::vector<double> x,
        std::vector<double> y,
        std::vector<double> weights = std::vector<double>())
{
    // sanity checks
    if (x.size() != y.size())
        throw std::runtime_error("x and y must have the same size.");
    if ((weights.size()) > 0 && (weights.size() != y.size()))
        throw std::runtime_error("x and y weights must have the same size.");

    // get inverse of permutation that brings x in ascending order
    std::vector<size_t> perm_x = utils::get_order(x);
    perm_x = utils::invert_permutation(perm_x);

    // sort x, y, and weights according to x, breaking ties with y
    utils::sort_all(x, y, weights);

    // get inverse of permutation that brings y in descending order
    std::vector<size_t> perm_y = utils::get_order(y, false);
    perm_y = utils::invert_permutation(perm_y);

    // sort y in descending order counting inversions
    std::vector<double> counts(y.size(), 0.0);
    merge_sort_count_per_element(y, weights, counts);

    // bring counts back in original order
    std::vector<double> counts_tmp = counts;
    for (size_t i = 0; i < counts.size(); i++)
        counts[i] = counts_tmp[perm_y[perm_x[i]]];

    return counts;
}

//! fast calculation of the weighted Hoeffdings's D.
//! @param x, y input data.
//! @param weights an optional vector of weights for the data.
inline double hoeffd(std::vector<double> x,
                     std::vector<double> y,
                     std::vector<double> weights = std::vector<double>())
{
    // 0. Check input sizes
    size_t n = x.size();
    if (y.size() != n)
        throw std::runtime_error("lengths of x and y must match.");
    bool weighted = (weights.size() > 0);
    if (weighted && (weights.size() != n))
        throw std::runtime_error("lengths of x, y, and weights must match.");

    if (weights.size() == 0)
        weights = std::vector<double>(n, 1);

    // 1. Compute (weighted) ranks
    std::vector<double> R_X = utils::compute_ranks(x, weights);
    std::vector<double> S_X = utils::compute_ranks(x, utils::pow(weights, 2));
    std::vector<double> R_Y = utils::compute_ranks(y, weights);
    std::vector<double> S_Y = utils::compute_ranks(y, utils::pow(weights, 2));

    // 2. Compute (weighted) bivariate ranks (number of points w/ both columns
    // less than the ith row).
    std::vector<double> R_XY(n), S_XY(n), T_XY(n), U_XY(n);
    R_XY = bivariate_rank(x, y, weights);
    S_XY = bivariate_rank(x, y, utils::pow(weights, 2));
    T_XY = bivariate_rank(x, y, utils::pow(weights, 3));
    U_XY = bivariate_rank(x, y, utils::pow(weights, 4));


    // 3. Compute (weighted) Hoeffdings' D
    double A_1 = 0.0, A_2 = 0.0, A_3 = 0.0;
    for (size_t i = 0; i < n; i++) {
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
    D += A_1 / (utils::perm_sum(weights, 3) * 6);
    D -= 2 * A_2 / (utils::perm_sum(weights, 4) * 24);
    D += A_3 / (utils::perm_sum(weights, 5) * 120);

    return 30.0 * D;
}

}
