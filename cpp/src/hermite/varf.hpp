#ifndef HERMITE_VARF_H
#define HERMITE_VARF_H

#include "hermite/types.hpp"
#include "boost/multi_array.hpp"

namespace hermite {

cube triple_products_1d(int degree);

template<typename T>
T varf(
        u_int degree,
        vec const & input,
        mat const & nodes,
        mat const & weights);

mat varfd(
        u_int dim,
        u_int degree,
        u_int direction,
        const mat & var);

spmat varfd(
        u_int dim,
        u_int degree,
        u_int direction,
        const spmat & var);

}

#endif
