#ifndef HERMITE_TRANSFORM_H
#define HERMITE_TRANSFORM_H

#include <string>
#include <vector>

#include "hermite/types.hpp"

namespace hermite {

std::vec discretize(
        std::s_func func,
        std::mat const & nodes,
        std::vec const & translation,
        std::mat const & dilation):

std::vec transform(
        std::u_int degree,
        std::vec const & f_grid,
        std::mat const & nodes,
        std::mat const & weights,
        bool forward);

double integrate(
        std::vec const & f_grid,
        std::mat const & nodes,
        std::mat const & weights);
}

#endif
