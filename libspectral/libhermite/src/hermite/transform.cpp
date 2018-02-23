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

#include "hermite/types.hpp"
#include "hermite/iterators.hpp"
#include "hermite/transform.hpp"

using namespace std;

namespace hermite {

void increment_multi_index(ivec & m, ivec const & upper_bounds)
{
    unsigned int dim = m.size();
    unsigned int i = dim - 1;
    while(m[i] == upper_bounds[i] - 1 && i>0) {
        m[i] = 0;
        i -= 1;
    }
    m[i] += 1;
}

vec hermite_expand(s_func func,
        unsigned int degree,
        cube const & nodes,
        cube const & weights,
        vec const & translation,
        mat const & dilation) {

    unsigned int dim = nodes[0].size();
    unsigned int n_products = nodes.size();

    using boost::math::binomial_coefficient;
    unsigned int n_polys = (unsigned int) binomial_coefficient<double> (degree + dim, dim);

    vec exp_coeffs(n_polys);

    unsigned int i,j,k,l;
    double* node = (double*) malloc(sizeof(double)*dim);
    double* mapped_node = (double*) malloc(sizeof(double)*dim);

    vec sq(degree + 1, 0);
    for (i = 0; i < sq.size(); i++)
    {
        sq[i] = sqrt(i);
    }

    for (i = 0; i < n_products; i++)
    {
        mat sub_nodes = nodes[i];
        mat sub_weights = weights[i];

        unsigned int n_points_tot = 1;
        ivec n_points(dim);
        for (j = 0; j < dim; j++)
        {
            n_points[j] = sub_nodes[j].size();
            n_points_tot *= n_points[j];
        }

        // Compute Hermite polynomials in each dimension
        cube herm_vals_1d(dim);
        for (j = 0; j < dim; ++j)
        {
            herm_vals_1d[j] = mat(n_points[j], vec(degree + 1, 0));
            for (k = 0; k < n_points[j]; k++)
            {
                double x = sub_nodes[j][k];
                herm_vals_1d[j][k][0] = 1;
                herm_vals_1d[j][k][1] = x;
                for (l = 1; l < degree; l++)
                {
                    herm_vals_1d[j][k][l+1] = (1/sq[l+1])*(x*herm_vals_1d[j][k][l]-sq[l]*herm_vals_1d[j][k][l-1]);
                }
            }
        }

        Hyper_cube_iterator iter_points(n_points);
        const ivec & p = iter_points.get();

        for (j = 0; j < n_points_tot; j++)
        {
            double weight = 1;
            for (k = 0; k < dim; k++)
            {
                node[k] = sub_nodes[k][p[k]];
                weight *= sub_weights[k][p[k]];
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

            Multi_index_iterator iter_degree(dim, degree);
            const ivec & m = iter_degree.get();

            for (k = 0; k < n_polys; k++)
            {
                double val_at_point = 1;
                for (l = 0; l < dim; l++)
                {
                    if (m[l] != 0)
                    {
                        val_at_point *=  herm_vals_1d[l][p[l]][k];
                    }
                }

                exp_coeffs[k] += func(mapped_node) * weight * val_at_point;
                iter_degree.increment();
            }
            iter_points.increment();
        }
    }

    return exp_coeffs;
}

double integrate_with_quad(s_func func,
        cube const & nodes,
        cube const & weights,
        vec const & translation,
        mat const & dilation) {

    vec integral = hermite_expand(func, 0, nodes, weights, translation, dilation);

    return integral[0];
}

void intern_function(string const & function_body) {

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


// ---- PYTHON WRAPPERS ----

double integrate_from_string(
        string const & function_body,
        cube const & nodes,
        cube const & weights,
        vec const & translation,
        mat const & dilation) {

    intern_function(function_body);
    string name = to_string(hash<string>()(function_body));
    string so_file = "/tmp/" + name + ".so";

    // Load function dynamically
    void *function_so = dlopen(so_file.c_str(), RTLD_NOW);
    s_func func = (s_func) dlsym(function_so, "toIntegrate");
    double result = integrate_with_quad(func, nodes, weights, translation, dilation);
    dlclose(function_so);
    return result;
}

vec hermite_expand_from_string(
        string const & function_body,
        int degree,
        cube const & nodes,
        cube const & weights,
        vec const & translation,
        mat const & dilation) {

    intern_function(function_body);
    string name = to_string(hash<string>()(function_body));
    string so_file = "/tmp/" + name + ".so";
    void *function_so = dlopen(so_file.c_str(), RTLD_NOW);
    s_func func = (s_func) dlsym(function_so, "toIntegrate");
    vec result = hermite_expand(func, degree, nodes, weights, translation, dilation);
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
    def("hermite_expand", hermite_expand_from_string);
}

}
