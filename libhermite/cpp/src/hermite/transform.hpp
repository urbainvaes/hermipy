#ifndef HERMITE_TRANSFORM_H
#define HERMITE_TRANSFORM_H

#include "hermite/types.hpp"

namespace hermite {

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
