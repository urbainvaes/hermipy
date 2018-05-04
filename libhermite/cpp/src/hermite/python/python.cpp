// #include <boost/python.hpp>
// #include <boost/python/numpy.hpp>
// #include <boost/core/ref.hpp>
// #include <boost/math/special_functions/binomial.hpp>
// #include <boost/python/suite/indexing/vector_indexing_suite.hpp>

// #include "hermite/discretize.hpp"
// #include "hermite/io.hpp"
// #include "hermite/iterators.hpp"
// #include "hermite/tensorize.hpp"
// #include "hermite/transform.hpp"
// #include "hermite/types.hpp"
// #include "hermite/varf.hpp"

// using namespace std;

// namespace hermite {

// namespace p = boost::python;
// namespace np = boost::python::numpy;

// c_mat contig_mat2(int rows, int cols)
// {
//     auto dims = boost::extents[rows][cols];
//     return boost::multi_array<double, 2>(dims);
// }


// mat test(int n)
// {
//     mat result(n, vec(n, 0.));
//     for (u_int i = 0; i < result.size(); i++)
//     {
//         for (u_int j = 0; j < result[0].size(); j++)
//         {
//             result[i][j] = (double) i - (double) j;
//         }
//     }
//     return result;
// }

// np::ndarray mat_to_numpy(mat const & input)
// {
//     u_int n_rows = input.size();
//     u_int n_cols = input[0].size();
//     p::tuple shape = p::make_tuple(n_rows, n_cols);
//     p::tuple stride = p::make_tuple(sizeof(double));
//     np::dtype dtype = np::dtype::get_builtin<double>();
//     p::object own;
//     np::ndarray converted = np::zeros(shape, dtype);

//     for (u_int i = 0; i < n_rows; i++)
//     {
//         shape = p::make_tuple(n_cols);
//         converted[i] = np::from_data(input[i].data(), dtype, shape, stride, own);
//     }
//     return converted;
// }

// mat to_cpp(np::ndarray input)
// {
//     u_int n_rows = input.size();
//     u_int n_cols = input[0].size();
//     p::tuple shape = p::make_tuple(n_rows, n_cols);
//     p::tuple stride = p::make_tuple(sizeof(double));
//     np::dtype dtype = np::dtype::get_builtin<double>();
//     p::object own;
//     np::ndarray converted = np::zeros(shape, dtype);

//     for (u_int i = 0; i < n_rows; i++)
//     {
//         shape = p::make_tuple(n_cols);
//         converted[i] = np::from_data(input[i].data(), dtype, shape, stride, own);
//     }
//     return converted;
// }

// np::ndarray cmat_to_numpy(c_mat const & input)
// {
//     u_int n_rows = input.shape()[0];
//     u_int n_cols = input.shape()[1];
//     p::tuple shape = p::make_tuple(n_rows, n_cols);
//     p::tuple strides = p::make_tuple(input.strides()[0]*sizeof(double),
//                                      input.strides()[1]*sizeof(double));
//     np::dtype dtype = np::dtype::get_builtin<double>();
//     p::object own;
//     np::ndarray converted = np::from_data(input.data(), dtype, shape, strides, own);
//     return converted;
// }

// c_mat test3(int n)
// {
//     c_mat result = contig_mat2(n, n);
//     for (u_int i = 0; i < result.shape()[0]; i++)
//     {
//         for (u_int j = 0; j < result.shape()[1]; j++)
//         {
//             result[i][j] = (double) i - (double) j;
//         }
//     }
//     return result;
// }

// np::ndarray test2(int n)
// {
//     c_mat result = test3(n);
//     return cmat_to_numpy(result);
// }


// // ---- PYTHON API ----


// }
