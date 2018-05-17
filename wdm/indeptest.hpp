// Copyright © 2018 Thomas Nagler
//
// This file is part of the wdm library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory
// or https://github.com/tnagler/wdmcpp/blob/master/LICENSE.

#pragma once

#include "utils.hpp"
#include "ktau.hpp"
#include "hoeffd.hpp"
#include "prho.hpp"
#include "srho.hpp"
#include "bbeta.hpp"
#include <boost/math/special_functions/atanh.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/math/interpolators/cubic_b_spline.hpp>
#include <boost/math/distributions/normal.hpp>
#include <random>
#include <limits>

namespace wdm {

//! calculates the (approximate) asymptotic distribution function of Hoeffding's
//! B (as in Blum, Kiefer, and Rosenblatt) under the null hypothesis of
//! independence.
//! @param B sample estimate of Hoeffding's B.
//! @param n the sample size.
inline double phoeffb(double B, size_t n) {
    B *= 0.5 * std::pow(boost::math::constants::pi<double>(), 4) * (n - 1);

    // obtain approximate p values by interpolation of tabulated values
    using namespace boost::math;
    double p;
    if ((B < 1.1) | (B > 8.6)) {
        p = std::min(1.0, std::exp(0.3885037 - 1.164879 * B));
        p = std::max(1e-12, p);
    } else if ((B >= 1.1) && (B < 5)) {
        std::vector<double> tab{
            0.5297, 0.4918, 0.4565, 0.4236, 0.3930, 0.3648, 0.3387, 0.3146,
            0.2924, 0.2719, 0.2530, 0.2355, 0.2194, 0.2045, 0.1908, 0.1781,
            0.1663, 0.1554, 0.1453, 0.1359, 0.1273, 0.1192, 0.1117, 0.1047,
            0.0982, 0.0921, 0.0864, 0.0812, 0.0762, 0.0716, 0.0673, 0.0633,
            0.0595, 0.0560, 0.0527, 0.0496, 0.0467, 0.0440, 0.0414, 0.0390,
            0.0368, 0.0347, 0.0327, 0.0308, 0.0291, 0.0274, 0.0259, 0.0244,
            0.0230, 0.0217, 0.0205, 0.0194, 0.0183, 0.0173, 0.0163, 0.0154,
            0.0145, 0.0137, 0.0130, 0.0123, 0.0116, 0.0110, 0.0104, 0.0098,
            0.0093, 0.0087, 0.0083, 0.0078, 0.0074, 0.0070, 0.0066, 0.0063,
            0.0059, 0.0056, 0.0053, 0.0050, 0.0047, 0.0045, 0.0042
        };
        cubic_b_spline<double> spline(tab.begin(), tab.end(), 1.1, 0.05);
        p = spline(B);
    } else {
        std::vector<double> tab{
            0.00025, 0.00014, 0.0008, 0.0005, 0.0003, 0.0002, 0.0001
        };
        cubic_b_spline<double> spline(tab.begin(), tab.end(), 5, 0.5);
        p = spline(B);
    }

    return p;
}

//! calculates the test statistic for indpendence tests
//! @param x, y input data.
//! @param method the dependence measure; possible values: `"prho"`, `"srho"`,
//!   `"ktau"`, `"hoeffd"`.
//! @param weights an optional vector of weights for the data.
inline double calculate_test_stat(
        const std::vector<double>& x,
        const std::vector<double>& y,
        std::string method,
        std::vector<double> weights = std::vector<double>())
{
    // determine effective sample size
    double n_eff;
    if (weights.size() == 0) {
        n_eff = static_cast<double>(x.size());
    } else {
        n_eff = wdm_utils::effective_sample_size(weights);
    }

    // calculate test statistic
    double stat;
    if (method == "hoeffd") {
        stat = hoeffd(x, y, weights) / 30.0 + 1.0 / (36.0 * n_eff);
    } else if (method == "ktau") {
        stat = ktau(x, y, weights);
        stat *= std::sqrt(9 * n_eff / 4);
    } else if (method == "prho") {
        stat = boost::math::atanh(prho(x, y, weights));
        stat *= std::sqrt(n_eff - 3);
    } else if (method == "srho") {
        stat = boost::math::atanh(srho(x, y, weights));
        stat *= std::sqrt((n_eff - 3) / 1.06);
    }  else if (method == "bbeta") {
        stat = bbeta(x, y, weights);
        stat *= std::sqrt(n_eff);
    } else {
        throw std::runtime_error("method not implemented.");
    }

    return stat;
}

//! calculates the asymptotic p-value.
//! @param stat value of the test statistic.
//! @param method the dependence measure; possible values: `"prho"`, `"srho"`,
//!   `"ktau"`, `"hoeffd"`.
//! @param n_eff effective sample size; only used for method `"hoeffd"`.
inline double calculate_asymptotic_p_val(double stat,
                                         std::string method,
                                         double n_eff = 0.0)
{
    double p_val;
    if (method == "hoeffd") {
        if (n_eff == 0.0)
            throw std::runtime_error("must provide n_eff for method 'hoeffd'.");
        p_val = phoeffb(stat, n_eff);
    } else {
        boost::math::normal norm_dist(0, 1);
        p_val = 2 * boost::math::cdf(norm_dist, -std::abs(stat));
    }

    return p_val;
}

//! calculates asymptotic p-values of independence tests based on (weighted)
//! dependence measures.
//! @param x, y input data.
//! @param method the dependence measure; possible values: `"prho"`, `"srho"`,
//!   `"ktau"`, `"bbeta"`, `"hoeffd"`.
//! @param weights an optional vector of weights for the data.
//! @return the p-value of the independence test.
inline double indeptest(
        const std::vector<double>& x,
        const std::vector<double>& y,
        std::string method,
        std::vector<double> weights = std::vector<double>())
{
    wdm_utils::check_sizes(x, y, weights);
    double stat = calculate_test_stat(x, y, method, weights);
    double n_eff = static_cast<double>(x.size());
    if (weights.size() > 0)
        n_eff = wdm_utils::effective_sample_size(weights);

    return calculate_asymptotic_p_val(stat, method, n_eff);
}

}
