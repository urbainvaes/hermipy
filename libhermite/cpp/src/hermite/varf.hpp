#ifndef HERMITE_VARF_H
#define HERMITE_VARF_H

#include "hermite/types.hpp"
#include "boost/multi_array.hpp"

namespace hermite {

typedef boost::multi_array<double, 2> c_mat;

std::cube triple_products_1d(int degree);

std::mat varf(
        std::u_int degree,
        std::vec const & input,
        std::mat const & nodes,
        std::mat const & weights);

std::mat varfd(
        std::u_int dim,
        std::u_int degree,
        std::u_int direction,
        std::mat const & var);
}

#endif
