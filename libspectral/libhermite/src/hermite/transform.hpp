#ifndef HERMITE_TRANSFORM_H
#define HERMITE_TRANSFORM_H

#include <string>
#include <vector>

#include "hermite/types.hpp"

namespace hermite {

std::vec hermite_expand(
        std::u_int degree,
        std::mat const & f_grid,
        std::cube const & nodes,
        std::cube const & weights);

double integrate_with_quad(
        std::mat const & f_grid,
        std::cube const & nodes,
        std::cube const & weights);

void intern_function(std::string const & function_body);

}

#endif
