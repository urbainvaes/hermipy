#include <dlfcn.h>
#include <fstream>

#include <boost/python.hpp>
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


// ---- PYTHON API ----
BOOST_PYTHON_MODULE(hermite_cpp)
{
    using namespace boost::python;

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
