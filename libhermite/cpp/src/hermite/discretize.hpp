#ifndef HERMITE_DISCRETIZE_H
#define HERMITE_DISCRETIZE_H

#include "hermite/types.hpp"

namespace hermite {

std::vec discretize(
        std::s_func func,
        std::mat const & nodes,
        std::vec const & translation,
        std::mat const & dilation);

void intern_function(std::string const & function_body);
std::vec discretize_from_string(
        std::string function_body,
        std::mat const & nodes,
        std::vec const & translation,
        std::mat const & dilation);

}
#endif
