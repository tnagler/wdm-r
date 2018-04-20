#pragma once

#include "utils.hpp"
#include "ktau.hpp"
#include "hoeffd.hpp"
#include "prho.hpp"
#include "srho.hpp"
#include <boost/math/special_functions/atanh.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/math/interpolators/cubic_b_spline.hpp>
#include <boost/math/distributions/normal.hpp>

namespace indep_test {

//! calculates the (approximate) asymptotic distribution function of Hoeffding's
//! D under the null hypothesis of independence using the method of Blum,
//! Kiefer, and Rosenblatt.
//! @param D sample estimate of Hoeffding's D.
//! @param n the sample size.
inline double phoeffd(double D, size_t n) {
    // transform to B statistic of Blum, Kiefer, and Rosenblatt
    D /= 30;
    D += 1.0 / (36.0 * n);
    D *= 0.5 * std::pow(boost::math::constants::pi<double>(), 4) * n;

    // obtain approximate p values by interpolation of tabulated values
    using namespace boost::math;
    double p;
    if ((D < 1.1) | (D > 8.6)) {
        p = std::min(1.0, std::exp(0.3885037 - 1.164879 * D));
        p = std::max(1e-12, p);
    } else if ((D >= 1.1) && (D < 5)) {
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
        p = spline(D);
    } else {
        std::vector<double> tab{
            0.00025, 0.00014, 0.0008, 0.0005, 0.0003, 0.0002, 0.0001
        };
        cubic_b_spline<double> spline(tab.begin(), tab.end(), 5, 0.5);
        p = spline(D);
    }

    return p;
}

//! calculates asymptotic p-values of independence tests based on (weighted)
//! dependence measures.
//! @param x, y input data.
//! @param weights an optional vector of weights for the data.
//! @param method the dependence measure; possible values: `"prho"`, `"srho"`,
//!   `"ktau"`, `"hoeffd"`.
//! @return the p-value of the independence test.
inline double indep_test_asymptotic(std::vector<double> x,
                                    std::vector<double> y,
                                    std::string method,
                                    std::vector<double> weights = std::vector<double>())
{
    utils::check_sizes(x, y, weights);
    double stat, p_val, n;

    // determine effective sample size
    if (weights.size() == 0)
        n = static_cast<double>(x.size());
    else
        n = utils::effective_sample_size(weights);

    // calculate test statistic and p-value
    if (method == "hoeffd") {
        stat = hoeffd::hoeffd(x, y, weights);
        p_val = phoeffd(stat, n);
    } else {
        if (method == "ktau") {
            stat = ktau::ktau(x, y, weights);
            stat *= std::sqrt(9 * n / 4);
        } else if (method == "prho") {
            stat = boost::math::atanh(prho::prho(x, y, weights));
            stat *= std::sqrt(n - 3);
        } else if (method == "srho") {
            stat = boost::math::atanh(srho::srho(x, y, weights));
            stat *= std::sqrt((n - 3) / 1.06);
        } else {
            throw std::runtime_error("method not impolemented.");
        }

        boost::math::normal normal_dist(0, 1);
        p_val = 2 * boost::math::cdf(normal_dist, -std::abs(stat));
    }

    return p_val;
}

}
