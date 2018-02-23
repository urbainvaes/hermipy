#ifndef HERMITE_TRANSFORM_H
#define HERMITE_TRANSFORM_H

#include <string>
#include <vector>

#include "hermite/types.hpp"

namespace hermite {

std::vec hermite_expand(
        std::u_int degree,
        std::vec const & f_grid,
        std::mat const & nodes,
        std::mat const & weights);

double integrate_with_quad(
        std::vec const & f_grid,
        std::mat const & nodes,
        std::mat const & weights);

void intern_function(std::string const & function_body);

}

#endif
