#pragma once

#include <algorithm>
#include <vector>
#include <numeric>

namespace utils {

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
    double n_eff = std::pow(utils::sum(weights), 2);
    n_eff /= utils::sum(utils::pow(weights, 2));

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
    std::sort(perm.begin(), perm.end(), utils::Sorter(x, ascending));

    return perm;
}


//! computes ranks (such that smallest element has rank 0), assigning average
//! ranks for ties.
//! @param x input vector.
//! @param ties_method `"min"` (default) assigns all tied values the minimum
//!   score; `"average"` assigns the average score.
//! @return a vector containing the ranks of each element in `x`.
std::vector<double> rank_scores(std::vector<double> x,
                                std::vector<double> weights = std::vector<double>(),
                                std::string ties_method = "min")
{
    size_t n = x.size();
    if (weights.size() == 0)
        weights = std::vector<double>(n, 1.0);

    std::vector<size_t> perm = get_order(x);

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

}
