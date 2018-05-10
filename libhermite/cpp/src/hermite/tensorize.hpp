#ifndef TENSORIZE_H
#define TENSORIZE_H

#include "hermite/types.hpp"
#include <boost/numeric/ublas/matrix_sparse.hpp>

namespace hermite {

std::vec tensorize_vecs(std::mat inputs);
std::mat tensorize_mats(std::cube inputs);
std::vec tensorize_vec(std::vec input, std::u_int dim, std::u_int dir);
std::mat tensorize_mat(std::mat input, std::u_int dim, std::u_int dir);
std::vec project_vec(std::vec input, std::u_int dim, std::u_int dir);
std::mat project_mat(std::mat input, std::u_int dim, std::u_int dir);

using namespace boost::numeric::ublas;
compressed_matrix<double, row_major> sp_tensorize_mats(std::vector<compressed_matrix<double, row_major>> inputs);

}

#endif

