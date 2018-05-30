#ifndef PROJECT_H
#define PROJECT_H

#include "hermite/types.hpp"

namespace hermite
{

vec project(const vec & input, u_int dim, u_int dir);
vec project(const vec & input, u_int dim, const ivec & dirs);

template <typename M> M project(const M & input, u_int dim, u_int dir);
template <typename M> M project(const M & input, u_int dim, const ivec & dirs);

}

#endif
