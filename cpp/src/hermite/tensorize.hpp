#ifndef TENSORIZE_H
#define TENSORIZE_H

#include "hermite/types.hpp"
#include <boost/numeric/ublas/matrix_sparse.hpp>

namespace hermite {

// Tensorize vector(s)
vec tensorize(const vec & input, u_int dim, u_int dir);
vec tensorize(const mat &);
vec tensorize(const mat & inputs, const imat & dirs);

// Tensorize matrix
template<typename T, typename M>
T tensorize(const M & input, u_int dim, u_int dir);

template <typename T, typename M>
T tensorize(const std::vector<M> & inputs);

template <typename T, typename S>
T tensorize(const std::vector<S> & inputs, const imat & dirs);

}

#endif
