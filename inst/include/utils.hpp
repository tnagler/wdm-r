#include <algorithm>
#include <vector>
#include <numeric>

namespace utils {


//! normalizes a vector of weights such that they sum to choose(n, 2), where
//! n is the total number of weights.
//! @param weights a vector of weights.
inline void normalize_weights(std::vector<double>& weights)
{
    size_t n = weights.size();
    if (n > 0) {
        double w_sum = 0.0, w2_sum = 0.0;
        for (size_t i = 0; i < weights.size(); i++) {
            w_sum  += weights[i];
            w2_sum += weights[i] * weights[i];
        }
        double norm = std::sqrt(0.5 * n * (n - 1));
        norm /= std::sqrt((w_sum * w_sum - w2_sum) / 2.0);
        for (size_t i = 0; i < weights.size(); i++) {
            weights[i] = weights[i] * norm;
        }
    }
}


//! inverts a permutation.
//! @param perm a permutation.
//! @return a vector containing the inverse permutation.
std::vector<size_t> invert_permuation(const std::vector<size_t>& perm)
{
    std::vector<size_t> inv_perm(perm.size());
    for (size_t i = 0; i < perm.size(); i++)
        inv_perm[perm[i]] = i;
    return inv_perm;
}

class Sorter {
public:
    Sorter(const std::vector<double>& x, bool ascending = true) : x_(x) {}

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

//! computes ranks.
//! @param x input vector.
//! @return a vector containing the ranks of each element in `x`.
std::vector<double> compute_ranks(std::vector<double> x,
                                  std::vector<double> weights = std::vector<double>())
{
    if (weights.size() == 0)
        weights = std::vector<double>(x.size(), 1.0);
    std::vector<size_t> perm = get_order(x);
    double w_acc = 0.0;
    for (size_t i = 0; i < x.size(); i++) {
        x[perm[i]] = w_acc;
        w_acc += weights[perm[i]];
    }

    return x;
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

inline std::vector<double> count_ties_per_element(
        const std::vector<double>& x,
        const std::vector<double>& weights)
{
    bool weighted = (weights.size() > 0);
    size_t n = x.size();
    std::vector<double> counts(n, 0.0);
    double w1 = 0.0, w2 = 0.0;
    size_t reps = 1;
    for (size_t i = 1; i < n; i++) {
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
                counts[i] += (w1 * w1 - w2) / 2.0;
            } else {
                counts[i] += reps * (reps - 1) / 2.0;
            }
            for (size_t k = 1; k <= reps; k++) {
                counts[i - k] = counts[i];
            }
            reps = 1;
        }
    }

    if (reps > 1) {
        if (weighted) {
            counts[n - 1] += (w1 * w1 - w2) / 2.0;
        } else {
            counts[n - 1] += reps * (reps - 1) / 2.0;
        }
        for (size_t k = 1; k <= reps; k++) {
            counts[n - 1 - k] = counts[n - 1];
        }
    }

    return counts;
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

inline std::vector<double> count_joint_ties_per_element(
        const std::vector<double>& x,
        const std::vector<double>& y,
        const std::vector<double>& weights)
{
    bool weighted = (weights.size() > 0);
    size_t n = x.size();
    std::vector<double> counts(n, 0.0);
    double w1 = 0.0, w2 = 0.0;
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
                counts[i] += (w1 * w1 - w2) / 2.0;
            } else {
                counts[i] += reps * (reps - 1) / 2.0;
            }
            for (size_t k = 1; k <= reps; k++)
                counts[i - k] = counts[i];
            reps = 1;
        }
    }

    if (reps > 1) {
        if (weighted) {
            counts[n - 1] += (w1 * w1 - w2) / 2.0;
        } else {
            counts[n - 1] += reps * (reps - 1) / 2.0;
        }
        for (size_t k = 1; k <= reps; k++)
            counts[n - 1 - k] = counts[n - 1];
    }

    return counts;
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
                counts[k] = counts2[j] + (w1_sum - w_acc) * weights2[j];
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
