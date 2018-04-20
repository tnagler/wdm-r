#pragma once

#include "utils.hpp"
#include "prho.hpp"

namespace srho {

//! fast calculation of the weighted Spearman's rho.
//! @param x, y input data.
//! @param weights an optional vector of weights for the data.
inline double srho(std::vector<double> x,
                   std::vector<double> y,
                   std::vector<double> weights = std::vector<double>())
{
    utils::check_sizes(x, y, weights);
    x = utils::rank_scores(x, weights);
    y = utils::rank_scores(y, weights);
    return prho::prho(x, y, weights);
}

}
