// Copyright © 2018 Thomas Nagler
//
// This file is part of the wdm library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory
// or https://github.com/tnagler/wdmcpp/blob/master/LICENSE.

#pragma once

#include <algorithm>
#include <vector>
#include <numeric>

namespace wdm_utils {

inline void check_sizes(const std::vector<double>& x,
                        const std::vector<double>& y,
                        const std::vector<double>& weights)
{
    if (y.size() != x.size())
        throw std::runtime_error("x and y must have the same size.");
    if ((weights.size() > 0) && (weights.size() != y.size()))
        throw std::runtime_error("x, y, and weights must have the same size.");
}


//! computes the nth power for all lements in a vector.
//! @param x the inpute vector.
//! @param n the exponent.
//! @return the vector x, but with all elements taken to the power n.
inline std::vector<double> pow(const std::vector<double>& x, size_t n)
{
    std::vector<double> res(x.size(), 1.0);
    if (n > 0) {
        for (size_t i = 0; i < x.size(); i++) {
            for (size_t j = 0; j < n; j++) {
                res[i] *= x[i];
            }
        }
    }

    return res;
}

//! sums all elements in a vector.
//! @param x the input vector.
inline double sum(const std::vector<double>& x) {
    double res = 0.0;
    for (size_t i = 0; i < x.size(); i++)
        res += x[i];
    return res;
}

//! computes the median of a vector.
//! @param x the input vector.
inline double median(std::vector<double> x) {
    size_t n = x.size();
    std::sort(x.begin(), x.end(), std::less<double>());
    double med;
    if (n % 2 == 0)
        med = 0.5 * (x[n / 2 - 1] + x[n / 2]);
    else
        med = x[n / 2];

    return med;
}

//! computes the sum of the products of all k-permutations of elements in a
//! vector using Newton's identities.
//! @param x the inpute vector.
//! @param k the order of the permutation.
inline double perm_sum(const std::vector<double>& x, size_t k) {
    if (k == 0)
        return 1.0;
    double s = 0;
    for (size_t i = 1; i <= k; i++)
        s += std::pow(-1, i - 1) * perm_sum(x, k - i) * sum(pow(x, i));
    return s / k;
}

//! computes the effective sample size from a sequence of weights.
//! @param the weight sequence.
inline double effective_sample_size(const std::vector<double>& weights) {
    double n_eff = std::pow(sum(weights), 2);
    n_eff /= sum(pow(weights, 2));

    return n_eff;
}

//! inverts a permutation.
//! @param perm a permutation.
//! @return a vector containing the inverse permutation.
std::vector<size_t> invert_permutation(const std::vector<size_t>& perm)
{
    std::vector<size_t> inv_perm(perm.size());
    for (size_t i = 0; i < perm.size(); i++)
        inv_perm[perm[i]] = i;
    return inv_perm;
}

//! Class providing an operator `(i, j)` that compares x[j] with x[j]; to be
//! used with `std::sort()`.
class Sorter {
public:
    //! Constructor
    //! @param x the reference vector.
    //! @param ascending whether to sort in ascending or descending order.
    Sorter(const std::vector<double>& x, bool ascending = true) :
    x_(x),
    ascending_(ascending) {}

    //! comparison operator.
    //! @param i first index.
    //! @param j second index.
    //! @return `(x_[i] < x_[j])` for ascending order; `(x_[i] > x_[j])` for
    //! descending order.
    inline bool operator()(size_t i, size_t j) const
    {
        if (ascending_)
            return (x_[i] < x_[j]);
        return (x_[i] > x_[j]);
    }

private:
    bool ascending_;
    const std::vector<double>& x_;
};


//! Class exporting an operator () comparing elements of x, breaking ties with y
class Sorter_with_tie_break {
public:
    Sorter_with_tie_break(const std::vector<double>& x,
                          const std::vector<double>& y)
        : x_(x), y_(y) {}

    inline bool operator()(size_t i, size_t j) const
    {
        if (x_[i] < x_[j]) {
            return true;
        } else if ((x_[i] == x_[j]) && (y_[i] < y_[j])) {
            return true;
        }
        return false;
    }
private:
    const std::vector<double>& x_;
    const std::vector<double>& y_;
} ;

//! computes the permutation that brings a vector into order.
//! @param x inpute vector.
//! @param ascending whether order ascendingly or descendingly.
std::vector<size_t> get_order(const std::vector<double>& x,
                              bool ascending = true)
{
    size_t n = x.size();
    std::vector<size_t> perm(n);
    for (size_t i = 0; i < n; i++)
        perm[i] = i;
    std::sort(perm.begin(), perm.end(), Sorter(x, ascending));

    return perm;
}

//! sorts x, y, and weights in x order; break ties in according to y.
//! @param x, y, weights input vectors.
inline void sort_all(std::vector<double>& x,
                     std::vector<double>& y,
                     std::vector<double>& weights)
{
    size_t n = x.size();
    std::vector<size_t> order(n);
    for (size_t i = 0; i < n; i++)
        order[i] = i;
    Sorter_with_tie_break sort_crit(x, y);
    std::sort(order.begin(), order.end(), sort_crit);

    std::vector<double> xx(n), yy(n);
    for (size_t i = 0; i < n; i++) {
        xx[i] = x[order[i]];
        yy[i] = y[order[i]];
    }

    // sort weights accordingly
    std::vector<double> w = weights;
    if (weights.size() > 0) {
        for (size_t i = 0; i < n; i++) {
            w[i] = weights[order[i]];
        }
    }

    x = xx;
    y = yy;
    weights = w;
}

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


//! sorts elements in a vector while counting inversions per element.
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

}
