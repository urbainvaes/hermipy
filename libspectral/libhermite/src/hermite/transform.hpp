#ifndef HERMITE_TRANSFORM_H
#define HERMITE_TRANSFORM_H

#include <string>
#include <vector>

#include "hermite/types.hpp"

namespace hermite {

void increment_multi_index(ivec & m, ivec const & upper_bounds);

vec hermite_expand(s_func func,
        unsigned int degree,
        cube const & nodes,
        cube const & weights,
        vec const & translation,
        mat const & dilation);

double integrate_with_quad(s_func func,
        cube const & nodes,
        cube const & weights,
        vec const & translation,
        mat const & dilation);

void intern_function(std::string const & function_body);

double integrate_from_string(
        std::string const & function_body,
        cube const & nodes,
        cube const & weights,
        vec const & translation,
        mat const & dilation);
}

#endif
