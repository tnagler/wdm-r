#pragma once

#include "utils.hpp"

namespace ktau {

//! count ties.
//! @param x a sorted input vector.
//! @param weights optionally, a vector of weights for the elements in `x`.
//! @param break_by second vector that was used to break ties (optional).
//! @return the number of (weighted) ties in `x`
inline double count_ties(const std::vector<double>& x,
                         const std::vector<double>& weights)
{
    bool weighted = (weights.size() > 0);
    double count = 0.0, w1 = 0.0, w2 = 0.0;
    size_t reps = 1;
    for (size_t i = 1; i < x.size(); i++) {
        if ((x[i] == x[i - 1])) {
            if (weighted) {
                if (reps == 1) {
                    w1 = weights[i - 1];
                    w2 = weights[i - 1] * weights[i - 1];
                }
                w1 += weights[i];
                w2 += weights[i] * weights[i];
            }
            reps++;
        } else if (reps > 1) {
            if (weighted) {
                count += (w1 * w1 - w2) / 2.0;
            } else {
                count += reps * (reps - 1) / 2.0;
            }
            reps = 1;
        }
    }

    if (reps > 1) {
        if (weighted) {
            count += (w1 * w1 - w2) / 2.0;
        } else {
            count += reps * (reps - 1) / 2.0;
        }
    }

    return count;
}

//! counts joint ties in two vectors.
//! @param x, y a input vectors that are sorted wrt `x` as first and `y` as
//!   secondary key.
//! @param weights optionally, a vector of weights for the elements in `x`.
//! @return the number of (weighted) joint ties in `x` and `y`.
inline double count_joint_ties(const std::vector<double>& x,
                               const std::vector<double>& y,
                               const std::vector<double>& weights)
{
    bool weighted = (weights.size() > 0);
    double count = 0.0, w1 = 0.0, w2 = 0.0;
    size_t reps = 1;
    size_t ref = 0;
    for (size_t i = 1; i < x.size(); i++) {
        if ((x[i] == x[i - 1]) && (y[i] == y[i - 1])) {
            if (weighted) {
                if (reps == 1) {
                    w1 = weights[i - 1];
                    w2 = weights[i - 1] * weights[i - 1];
                }
                w1 += weights[i];
                w2 += weights[i] * weights[i];
            }
            reps++;
        } else if (reps > 1) {
            if (weighted) {
                count += (w1 * w1 - w2) / 2.0;
            } else {
                count += reps * (reps - 1) / 2.0;
            }
            reps = 1;
        }
    }

    if (reps > 1) {
        if (weighted) {
            count += (w1 * w1 - w2) / 2.0;
        } else {
            count += reps * (reps - 1) / 2.0;
        }
    }

    return count;
}

//! merge sort for a pair of vectors, counting inversions.
//! @param vec container for the sorted elements.
//! @param vec1, vec2 sorted input vectors to be merged.
//! @param weights container for the weights corresponding to sorted elements
//!   in `vec`; can be empty for unweighted counts.
//! @param weights1, weights2 weights corresponding to input vectors `vec1`,
//!   `vec2`; can be empty for unweighted counts.
//! @param count counter to which the (weighted) number of inversions is added.
inline void merge(std::vector<double>& vec,
                  const std::vector<double>& vec1,
                  const std::vector<double>& vec2,
                  std::vector<double>& weights,
                  const std::vector<double>& weights1,
                  const std::vector<double>& weights2,
                  double& count)
{
    double w_acc = 0.0, w1_sum = 0.0;
    bool weighted = (weights.size() > 0);
    if (weighted) {
        for (size_t i = 0; i < weights1.size(); i++)
            w1_sum += weights1[i];
    }
    size_t i, j, k;
    for (i = 0, j = 0, k = 0; i < vec1.size() && j < vec2.size(); k++) {
        if (vec1[i] <= vec2[j]) {
            vec[k] = vec1[i];
            if (weighted) {
                weights[k] = weights1[i];
                w_acc += weights1[i];
            }
            i++;
        } else {
            vec[k] = vec2[j];
            if (weighted) {
                weights[k] = weights2[j];
                count += weights2[j] * (w1_sum - w_acc);
            } else {
                count += vec1.size() - i;
            }
            j++;
        }
    }

    while (i < vec1.size()) {
        vec[k] = vec1[i];
        if (weighted)
            weights[k] = weights1[i];
        k++;
        i++;
    }

    while (j < vec2.size()) {
        vec[k] = vec2[j];
        if (weighted)
            weights[k] = weights2[j];
        k++;
        j++;
    }
}

//! sorting elements in a vector while counting inversions.
//! @param vec the vector to be sorted.
//! @param weights vector of weights corresponding to `vec`; can be empty for
//!   unweighted counts.
//! @param count counter to which the (weighted) number of inversions are added.
inline void merge_sort(std::vector<double>& vec,
                       std::vector<double>& weights,
                       double& count)
{
    if (vec.size() > 1) {
        size_t n = vec.size();
        std::vector<double> vec1(vec.begin(), vec.begin() + n / 2);
        std::vector<double> vec2(vec.begin() + n / 2, vec.end());

        n = weights.size();
        std::vector<double> weights1(weights.begin(), weights.begin() + n / 2);
        std::vector<double> weights2(weights.begin() + n / 2, weights.end());

        merge_sort(vec1, weights1, count);
        merge_sort(vec2, weights2, count);
        merge(vec, vec1, vec2, weights, weights1, weights2, count);
    }
}

//! fast calculation of the weighted Kendall's tau.
//! @param x, y input data.
//! @param weights an optional vector of weights for the data.
inline double ktau(std::vector<double> x,
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

    // 1.1 Sort x, y, and weights in x order; break ties in according to y.
    utils::sort_all(x, y, weights);

    // 1.2 Count pairs of tied x and simultaneous ties in x and y.
    double ties_x = count_ties(x, weights);
    double ties_both = count_joint_ties(x, y, weights);

    // 2.1 Sort y again and count exchanges (= number of discordant pairs).
    double num_d = 0.0;
    merge_sort(y, weights, num_d);

    // 2.2 Count pairs of tied y.
    double ties_y = count_ties(y, weights);

    // 3. Calculate Kendall's tau.
    double num_pairs = utils::perm_sum(weights, 2);
    double num_c = num_pairs - (num_d + ties_x + ties_y - ties_both);
    double tau = num_c - num_d;
    tau /= std::sqrt((num_pairs - ties_x) * (num_pairs - ties_y));

    return tau;
}

}
