#ifndef TENSORIZE_H
#define TENSORIZE_H

#include "hermite/types.hpp"
#include <boost/numeric/ublas/matrix_sparse.hpp>

namespace hermite {

// Tensorize vector(s)
std::vec tensorize(const std::mat &);
std::vec tensorize(const std::vec & input, std::u_int dim, std::u_int dir);

// Project vector or matrix
std::vec project(const std::vec & input, std::u_int dim, std::u_int dir);
std::mat project(const std::mat & input, std::u_int dim, std::u_int dir);

// Tensorize matrix
template<typename T> T tensorize(const std::mat & input, std::u_int dim, std::u_int dir);
template<typename T> T tensorize(const std::cube & inputs);

}

#endif

