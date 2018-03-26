#ifndef HERMITE_VARF_H
#define HERMITE_VARF_H

#include "hermite/types.hpp"

namespace hermite {

std::cube triple_products_1d(int degree);

std::mat varf(
        std::u_int degree,
        std::vec const & input,
        std::mat const & nodes,
        std::mat const & weights);
}
#endif
