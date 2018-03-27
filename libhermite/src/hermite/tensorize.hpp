#ifndef TENSORIZE_H
#define TENSORIZE_H

#include "hermite/types.hpp"

namespace hermite {

std::vec tensorize_vec(std::vec input, std::u_int dim, std::u_int dir);
std::mat tensorize_mat(std::mat input, std::u_int dim, std::u_int dir);

}

#endif

