#include <cmath>
#include <dlfcn.h>
#include <functional>
#include <numeric>
#include <iostream>
#include <fstream>
#include <tuple>
#include <utility>
#include <vector>

#include <boost/function.hpp>
#include <boost/python.hpp>
#include <boost/core/ref.hpp>
#include <boost/math/special_functions/binomial.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>


// #include "hermite/helpers/templates.hpp"
// #include "hermite/helpers/combinatorics.hpp"
// #include "hermite/integrator.hpp"
// #include "hermite/quadrature.hpp"

using namespace std;

namespace hermite {

// Arrays
typedef std::vector<double> vec;
typedef std::vector<vec> mat;
typedef std::vector<mat> cube;

typedef std::vector<int> ivec;

// Functions
// typedef boost::function<double(vec const &)> s_func;
typedef double (*s_func)(double*);

void increment_multi_index(ivec & m, ivec const & upper_bounds)
{
    unsigned int dim = m.size();
    unsigned int i = dim - 1;
    while(m[i] == upper_bounds[i] - 1 && i>0) {
        i -= 1;
    }
    m[i] += 1;
    for (unsigned int j = i + 1; j < dim; j++)
    {
        m[j] = 0;
    }
}

// vec hermite_expand(s_func func,

double integrate_with_quad(s_func func,
        cube const & nodes,
        cube const & weights,
        vec const & translation,
        mat const & dilation) {

    unsigned int dim = nodes[0].size();
    unsigned int n_products = nodes.size();

    double result = 0.;

    unsigned int i,j,k,l;
    double* node = (double*) malloc(sizeof(double)*3);
    double* mapped_node = (double*) malloc(sizeof(double)*3);
    // vec node(dim);
    // vec mapped_node(dim);

    for (i = 0; i < n_products; ++i)
    {
        mat sub_nodes = nodes[i];
        mat sub_weights = weights[i];
        ivec n_points(dim);
        unsigned int n_points_tot = 1;

        for (j = 0; j < dim; j++)
        {
            n_points[j] = sub_nodes[j].size();
            n_points_tot *= n_points[j];
        }

        ivec m(dim, 0);
        for (j = 0; j < n_points_tot; j++)
        {
            double weight = 1;
            for (k = 0; k < dim; k++)
            {
                node[k] = sub_nodes[k][m[k]];
                weight *= sub_weights[k][m[k]];
            }

            for (k = 0; k < dim; k++)
            {
                mapped_node[k] = 0;
                for (l = 0; l < dim; l++)
                {
                    mapped_node[k] += dilation[k][l]*node[l];
                }
                mapped_node[k] += translation[k];
            }

            result += weight * func(mapped_node);
            increment_multi_index(m, n_points);
        }
    }
    return result;
}

double integrate_from_string(
        string const & function_body,
        cube const & nodes,
        cube const & weights,
        vec const & translation,
        mat const & dilation) {

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

    // Load function dynamically
    void *function_so = dlopen(so_file.c_str(), RTLD_NOW);
    s_func func = (s_func) dlsym(function_so, "toIntegrate");
    double result = integrate_with_quad(func, nodes, weights, translation, dilation);
    dlclose(function_so);
    return result;
}


// ---- PYTHON API ----
BOOST_PYTHON_MODULE(hermite)
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

    def("integrate_from_string", integrate_from_string);
}

}
