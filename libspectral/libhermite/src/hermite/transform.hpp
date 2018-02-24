#ifndef HERMITE_TRANSFORM_H
#define HERMITE_TRANSFORM_H

#include <string>
#include <vector>

#include "hermite/types.hpp"

namespace hermite {

std::vec hermite_transform(
        std::u_int degree,
        std::vec const & f_grid,
        std::mat const & nodes,
        std::mat const & weights,
        bool forward);

double integrate_with_quad(
        std::vec const & f_grid,
        std::mat const & nodes,
        std::mat const & weights);
}

#endif
