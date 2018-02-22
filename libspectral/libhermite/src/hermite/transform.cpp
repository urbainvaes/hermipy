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
typedef std::vector<ivec> imat;
typedef std::vector<imat> icube;

// Functions
typedef std::function<double(vec const &)> s_func;

void increment_multi_index(ivec &m, ivec const & upper_bounds)
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

double integrate_with_quad(s_func const & func, cube const & nodes, cube const & weights) {

    unsigned int i,j,k;
    unsigned int dim = nodes[0].size();
    unsigned int n_products = nodes.size();

    double result = 0.;

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
            vec node(dim);
            double weight = 1;
            for (k = 0; k < dim; k++)
            {
                node[k] = sub_nodes[k][m[k]];
                weight *= sub_weights[k][m[k]];
            }
            result += weight * func(node);
            increment_multi_index(m, n_points);
        }
    }
    return result;
}

double integrate_from_string(string const& function_body, cube nodes, cube weights) {

    // Write function to file
    ofstream helper_file;
    helper_file.open("/tmp/helper_function.cpp");
    helper_file << "#include <vector>\n#include <cmath>\n";
    helper_file << "extern \"C\" double toIntegrate(std::vector<double> v) {\n";
    helper_file << "    return " << function_body << ";\n}";
    helper_file.close();

    // Compile file
    system("c++ /tmp/helper_function.cpp -o /tmp/helper_function.so -shared -fPIC");

    // Load function dynamically
    typedef double (*vec_func)(vec);
    void *function_so = dlopen("/tmp/helper_function.so", RTLD_NOW);
    vec_func func = (vec_func) dlsym(function_so, "toIntegrate");
    double result = integrate_with_quad(func, nodes, weights);
    dlclose(function_so);
    return result;
}

// // ---- PYTHON WRAPPERS ----
// double Quad::integrate_wrapper(boost::python::object const& func) {
//     std::function<double(vec const&)> lambda;
//     switch (this->nodes[0].size()) {
//         case 1: lambda = [&func](vec const& v) -> double { return boost::python::extract<double>(func (v[0])); }; break;
//         case 2: lambda = [&func](vec const& v) -> double { return boost::python::extract<double>(func(v[0], v[1])); }; break;
//         case 3: lambda = [&func](vec const& v) -> double { return boost::python::extract<double>(func(v[0], v[1], v[2])); }; break;
//         default: cout << "Dimension must be 1, 2, or 3" << endl; exit(0);
//     }
//     return integrate(lambda);
// }

// ---- PYTHON API ----
BOOST_PYTHON_MODULE(hermite)
{
    using namespace boost::python;

    class_<vec>("cpp_vector")
        .def(vector_indexing_suite<vec>())
        ;

    class_<mat>("cpp_mat")
        .def(vector_indexing_suite<mat>())
        ;

    class_<cube>("cpp_cube")
        .def(vector_indexing_suite<cube>())
        ;


    def("integrate_from_string", integrate_from_string);

    // class_<Quad>("Quad", init<int,int>())
    //     .def("integrate", &Quad::integrate_wrapper)
    //     .def("integrate_from_string", &integrate_from_string)
    //     .def_readonly("nodes", &Quad::nodes)
    //     .def_readonly("weights", &Quad::weights)
    //     ;
}

}
