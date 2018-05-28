#ifndef SPARSE_H
#define SPARSE_H

#include <boost/numeric/ublas/matrix_sparse.hpp>
#include "hermite/types.hpp"

namespace bnu = boost::numeric::ublas;
namespace hermite {

mat row_col_val(const spmat & input);
mat full(const spmat & input);
spmat to_spmat(const mat & input);

spmat to_spmat(
        const vec & data,
        const ivec & indices,
        const ivec & indptr,
        u_int size1,
        u_int size2);

}

#endif
