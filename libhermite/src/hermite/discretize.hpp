#ifndef HERMITE_DISCRETIZE_H
#define HERMITE_DISCRETIZE_H

#include "hermite/types.hpp"

namespace hermite {

std::vec discretize(
        std::s_func func,
        std::mat const & nodes,
        std::vec const & translation,
        std::mat const & dilation);

}
#endif
