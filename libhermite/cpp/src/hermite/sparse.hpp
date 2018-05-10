#ifndef SPARSE_H
#define SPARSE_H

#include <boost/numeric/ublas/matrix_sparse.hpp>
#include "hermite/types.hpp"

namespace bnu = boost::numeric::ublas;
namespace hermite {

std::mat spmat_to_mat(bnu::compressed_matrix<double, bnu::row_major> input);
bnu::compressed_matrix<double, bnu::row_major> mat_to_spmat(std::mat input);

}

#endif
