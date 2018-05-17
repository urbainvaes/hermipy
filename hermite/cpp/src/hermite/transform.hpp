#ifndef HERMITE_TRANSFORM_H
#define HERMITE_TRANSFORM_H

#include "hermite/types.hpp"

namespace hermite {

vec transform(
        u_int degree,
        vec const & f_grid,
        mat const & nodes,
        mat const & weights,
        bool forward);

double integrate(
        vec const & f_grid,
        mat const & nodes,
        mat const & weights);
}

#endif
