#ifndef HERMITE_DISCRETIZE_H
#define HERMITE_DISCRETIZE_H

#include "hermite/types.hpp"

namespace hermite {

vec discretize(
        s_func func,
        mat const & nodes,
        vec const & translation,
        mat const & dilation);

void intern_function(std::string const & function_body);
vec discretize_from_string(
        std::string function_body,
        mat const & nodes,
        vec const & translation,
        mat const & dilation);

}
#endif
