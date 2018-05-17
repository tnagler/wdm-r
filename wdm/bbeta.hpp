#pragma once

#include "utils.hpp"

namespace wdm {

//! calculates the weighted Blomqvists's beta.
//! @param x, y input data.
//! @param weights an optional vector of weights for the data.
inline double bbeta(const std::vector<double>& x,
                    const std::vector<double>& y,
                    std::vector<double> weights = std::vector<double>())
{
    wdm_utils::check_sizes(x, y, weights);
    size_t n = x.size();
    if (weights.size() == 0)
        weights = std::vector<double>(n, 1.0);

    // find the medians
    double med_x = wdm_utils::median(x);
    double med_y = wdm_utils::median(y);

    // count elements in lower left and upper right quadrants
    double w_acc{0.0};
    for (size_t i = 0; i < n; i++) {
        if ((x[i] <= med_x) && (y[i] <= med_y))
            w_acc += weights[i];
        else if ((x[i] > med_x) && (y[i] > med_y))
            w_acc += weights[i];
    }

    return 2 * w_acc / wdm_utils::sum(weights) - 1;
}

}
