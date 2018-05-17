#pragma once

#include "utils.hpp"
#include "ranks.hpp"
#include "prho.hpp"

namespace wdm {

//! fast calculation of the weighted Spearman's rho.
//! @param x, y input data.
//! @param weights an optional vector of weights for the data.
inline double srho(std::vector<double> x,
                   std::vector<double> y,
                   std::vector<double> weights = std::vector<double>())
{
    wdm_utils::check_sizes(x, y, weights);
    x = rank_scores(x, weights, "average");
    y = rank_scores(y, weights, "average");
    return prho(x, y, weights);
}

}
