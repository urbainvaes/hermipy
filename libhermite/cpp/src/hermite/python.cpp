#include <dlfcn.h>
#include <fstream>

#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include <boost/core/ref.hpp>
#include <boost/math/special_functions/binomial.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

#include "hermite/discretize.hpp"
#include "hermite/io.hpp"
#include "hermite/iterators.hpp"
#include "hermite/tensorize.hpp"
#include "hermite/transform.hpp"
#include "hermite/types.hpp"
#include "hermite/varf.hpp"

using namespace std;

namespace hermite {

void intern_function(string const & function_body)
{
    string name = to_string(hash<string>()(function_body));
    string cpp_file = "/tmp/" + name + ".cpp";
    string so_file = "/tmp/" + name + ".so";
    ifstream test_exists(so_file.c_str());

    if(! test_exists.good()) {
         ofstream write_function;
         write_function.open(cpp_file);
         write_function << "#include <vector>\n#include <cmath>\n";
         write_function << "extern \"C\" double toIntegrate(double *v) {\n";
         write_function << "    return " << function_body << ";\n}";
         write_function.close();

        // Compile file
        string command = "c++ " + cpp_file + " -o " + so_file + " -O3 -Ofast -shared -fPIC";
        system(command.c_str());
    }
}

vec discretize_from_string(
        string function_body,
        mat const & nodes,
        vec const & translation,
        mat const & dilation)
{
    intern_function(function_body);
    string name = to_string(hash<string>()(function_body));
    string so_file = "/tmp/" + name + ".so";
    void *function_so = dlopen(so_file.c_str(), RTLD_NOW);
    s_func func = (s_func) dlsym(function_so, "toIntegrate");
    vec result = discretize(func, nodes, translation, dilation);
    dlclose(function_so);
    return result;
}

namespace p = boost::python;
namespace np = boost::python::numpy;

c_mat contig_mat2(int rows, int cols)
{
    auto dims = boost::extents[rows][cols];
    return boost::multi_array<double, 2>(dims);
}


mat test(int n)
{
    mat result(n, vec(n, 0.));
    for (u_int i = 0; i < result.size(); i++)
    {
        for (u_int j = 0; j < result[0].size(); j++)
        {
            result[i][j] = (double) i - (double) j;
        }
    }
    return result;
}

np::ndarray mat_to_numpy(mat const & input)
{
    u_int n_rows = input.size();
    u_int n_cols = input[0].size();
    p::tuple shape = p::make_tuple(n_rows, n_cols);
    p::tuple stride = p::make_tuple(sizeof(double));
    np::dtype dtype = np::dtype::get_builtin<double>();
    p::object own;
    np::ndarray converted = np::zeros(shape, dtype);

    for (u_int i = 0; i < n_rows; i++)
    {
        shape = p::make_tuple(n_cols);
        converted[i] = np::from_data(input[i].data(), dtype, shape, stride, own);
    }
    return converted;
}

np::ndarray cmat_to_numpy(c_mat const & input)
{
    u_int n_rows = input.shape()[0];
    u_int n_cols = input.shape()[1];
    p::tuple shape = p::make_tuple(n_rows, n_cols);
    p::tuple strides = p::make_tuple(input.strides()[0]*sizeof(double),
                                     input.strides()[1]*sizeof(double));
    np::dtype dtype = np::dtype::get_builtin<double>();
    p::object own;
    np::ndarray converted = np::from_data(input.data(), dtype, shape, strides, own);
    return converted;
}

c_mat test3(int n)
{
    c_mat result = contig_mat2(n, n);
    for (u_int i = 0; i < result.shape()[0]; i++)
    {
        for (u_int j = 0; j < result.shape()[1]; j++)
        {
            result[i][j] = (double) i - (double) j;
        }
    }
    return result;
}

np::ndarray test2(int n)
{
    c_mat result = test3(n);
    return cmat_to_numpy(result);
}


// ---- PYTHON API ----
BOOST_PYTHON_MODULE(hermite_cpp)
{
    using namespace boost::python;

    // Initialize numpy
    Py_Initialize();
    boost::python::numpy::initialize();

    class_<vec>("double_vec")
        .def(vector_indexing_suite<vec>())
        ;

    class_<mat>("double_mat")
        .def(vector_indexing_suite<mat>())
        ;

    class_<cube>("double_cube")
        .def(vector_indexing_suite<cube>())
        ;

    class_<ivec>("int_vec")
        .def(vector_indexing_suite<ivec>())
        ;

    class_<imat>("int_mat")
        .def(vector_indexing_suite<imat>())
        ;

    class_<c_mat>("contiguous_mat");

    def("test", test);
    def("test2", test2);
    def("test3", test3);
    def("to_numpy", mat_to_numpy);
    def("to_numpy", cmat_to_numpy);
    def("discretize", discretize_from_string);
    def("integrate", integrate);
    def("list_cube_indices", list_cube_indices);
    def("list_multi_indices", list_multi_indices);
    def("project", project_vec);
    def("project", project_mat);
    def("tensorize", tensorize_vecs);
    def("tensorize", tensorize_mats);
    def("tensorize", tensorize_vec);
    def("tensorize", tensorize_mat);
    def("transform", transform);
    def("triple_products", triple_products_1d);
    def("varf", varf);
    def("varfd", varfd);
}

}
