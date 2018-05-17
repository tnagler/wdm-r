#pragma once

#include "wdm/ktau.hpp"
#include "wdm/hoeffd.hpp"
#include "wdm/prho.hpp"
#include "wdm/srho.hpp"
#include "wdm/bbeta.hpp"
#include "wdm/indeptest.hpp"

//! Weighted dependence measures
namespace wdm {

//! calculates (weighted) dependence measures.
//! @param x, y input data.
//! @param method the dependence measure; possible values: `"prho"`, `"srho"`,
//!   `"ktau"`, `"bbeta"`, `"hoeffd"`.
//! @param weights an optional vector of weights for the data.
//! @return the dependence measure
double wdm(const std::vector<double>& x,
           const std::vector<double>& y,
           std::string method,
           std::vector<double> weights = std::vector<double>())
{
    if (method == "hoeffd") {
        return hoeffd(x, y, weights);
    } else if (method == "ktau") {
        return ktau(x, y, weights);
    } else if (method == "prho") {
        return prho(x, y, weights);
    } else if (method == "srho") {
        return srho(x, y, weights);
    }  else if (method == "bbeta") {
        return bbeta(x, y, weights);
    } else {
        throw std::runtime_error("method not implemented.");
    }
}

}
