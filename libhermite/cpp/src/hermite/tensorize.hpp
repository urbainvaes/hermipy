#ifndef TENSORIZE_H
#define TENSORIZE_H

#include "hermite/types.hpp"
#include <boost/numeric/ublas/matrix_sparse.hpp>

namespace hermite {

// Tensorize vector(s)
vec tensorize(const mat &);
vec tensorize(const vec & input, u_int dim, u_int dir);

// Project vector or matrix
vec project(const vec & input, u_int dim, u_int dir);
mat project(const mat & input, u_int dim, u_int dir);

// Tensorize matrix
template<typename T, typename M> 
T tensorize(const M & input, u_int dim, u_int dir);

template <typename T, typename M>
T tensorize(const std::vector<M> & inputs);

}

#endif

