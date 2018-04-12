#include "utils.hpp"

namespace fasthoeffd {

inline double perm_sum(const std::vector<double>& x, size_t k) {
    if (k == 0)
        return 1.0;
    double s = 0;
    for (size_t i = 1; i <= k; i++) {
        double xi_sum = 0.0;
        for (size_t j = 0; j < x.size(); j++)
            xi_sum += std::pow(x[j], i);
        s += std::pow(-1, i - 1) * perm_sum(x, k - i) * xi_sum;
    }
    return s / k;
}


inline double fast_hoeffd(std::vector<double> x,
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

    // 1. Compute the ranks
    std::vector<double> r0 = utils::compute_ranks(x, weights);
    std::vector<double> r1 = utils::compute_ranks(y, weights);
    // for (size_t i = 0; i < n; i++)
    //     std::cout << r0[i] << ", " << r1[i] << std::endl;

    // 2. Sort all vectors according to first variable, breaking ties with second.
    utils::sort_all(r0, r1, weights);

    // 3. Compute Qi (number of points w/ both columns less than the ith row).
    // 3.1 Find permutation corresponding to a sorting in descending order
    std::vector<size_t> perm = utils::get_order(r1, false);

    // 3.2 Q is computetd by sorting the second variable again in descending
    // order while counting inversions.
    std::vector<double> Q(n, 0.0);
    auto ones = std::vector<double>(n, 1);
    if (weights.size() == 0)
        weights = ones;
    utils::merge_sort_count_per_element(r1, weights, Q);

    // 4. Compute Hoeffdings' D
    double A = 0.0, A1 = 0.0, A2 = 0.0, A3 = 0.0,
        B = 0.0, B1 = 0.0, B2 = 0.0, B3 = 0.0,
        C = 0.0, CC = 0.0;
    for (size_t i = 0; i < n; i++) {
        A += (r0[perm[i]] - 1) * (r0[perm[i]] - 1) * (r1[i]) * (r1[i]) * weights[i];
        A1 -= (r0[perm[i]]) * (r1[i]) * (r1[i]) * weights[i];
        A2 -= (r0[perm[i]]) * (r0[perm[i]]) * (r1[i]) * weights[i];
        A3 += (r0[perm[i]]) * (r1[i]) * weights[i];
        B += r0[perm[i]] * r1[i] * Q[i] * weights[i];
        B1 -= r1[i] * Q[i] * weights[i];
        B2 -= r0[perm[i]] * Q[i] * weights[i];
        B3 += Q[i] * weights[i];
        C += Q[i] * Q[i] * weights[i];
        CC -= Q[i] * weights[i];
    }

    for (size_t i = 0; i < n; i++) {
        A += r0[perm[i]] * r0[perm[i]] * r1[i] * r1[i] * weights[i];
        B += r0[perm[i]] * r1[i] * Q[i] * weights[i];
        C += Q[i] * Q[i] * weights[i];
    }

    return (A - 2 * B + C) / perm_sum(weights, 5) / 120;

    double D = A - 2 * (n - 2) * B + (n - 2) * (n - 3) * C;
    D /= perm_sum(weights, 5) * 120;


    return 30.0 * D;
}

}