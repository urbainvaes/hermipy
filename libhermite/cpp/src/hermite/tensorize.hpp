#ifndef TENSORIZE_H
#define TENSORIZE_H

#include "hermite/types.hpp"
#include <boost/numeric/ublas/matrix_sparse.hpp>

namespace hermite {

std::vec tensorize(const std::mat &);
std::mat tensorize(const std::cube &);
std::vec tensorize(const std::vec & input, std::u_int dim, std::u_int dir);
std::mat tensorize(const std::mat & input, std::u_int dim, std::u_int dir);
std::vec project(const std::vec & input, std::u_int dim, std::u_int dir);
std::mat project(const std::mat & input, std::u_int dim, std::u_int dir);

// using namespace boost::numeric::ublas;
boost::spmat tensorize(const std::vector<boost::spmat> & inputs);

}

#endif

