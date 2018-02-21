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


#include "hermite/helpers/templates.hpp"
#include "hermite/helpers/combinatorics.hpp"
#include "hermite/integrator.hpp"
#include "hermite/quadrature.hpp"

typedef std::vector< std::vector< std::vector<double> > > cube;
typedef std::vector< std::vector<double> > mat;
typedef std::vector<double> vec;

using namespace std;

namespace hermite {


    // Compute the nodes and weights resulting from the product of quadratures
    void quad_prod(vector<int> sizes, mat& nodes, vec& weights) {

        // Dimension
        int dim = sizes.size();

        // Number of nodes of product
        int n = 1;
        for (int i = 0; i < dim; ++i)
            n *= sizes[i];

        // Store nodes and weights of all the quadrature rules used
        mat all_nodes(dim);
        mat all_weights(dim);

        for (int i = 0; i < dim; ++i)
            get_gh_quadrature(sizes[i], all_nodes[i], all_weights[i]);

        // Initialize
        nodes = mat (n, vec(dim, 0.));
        weights = vec (n, 1.);

        // Compute nodes and weights of product
        for (int i = 0; i < n; ++i) {
            for (int j = 0, aux = i; j < dim; ++j) {
                nodes[i][j] = all_nodes[j][aux % sizes[j]];
                weights[i] *= all_weights[j][aux % sizes[j]];
                aux = aux / sizes[j];
            }
        }
    }

    double Quad::integrate_from_string(string const& function_body) {

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
        double result = integrate(func);
        dlclose(function_so);
        return result;
    }

    // ---- PYTHON WRAPPERS ----
    double Quad::integrate_wrapper(boost::python::object const& func) {
        std::function<double(vec const&)> lambda;
        switch (this->nodes[0].size()) {
            case 1: lambda = [&func](vec const& v) -> double { return boost::python::extract<double>(func (v[0])); }; break;
            case 2: lambda = [&func](vec const& v) -> double { return boost::python::extract<double>(func(v[0], v[1])); }; break;
            case 3: lambda = [&func](vec const& v) -> double { return boost::python::extract<double>(func(v[0], v[1], v[2])); }; break;
            default: cout << "Dimension must be 1, 2, or 3" << endl; exit(0);
        }
        return integrate(lambda);
    }

    // ---- EXPOSE TO PYTHON ----
    BOOST_PYTHON_MODULE(hermite)
    {
        using namespace boost::python;

        class_<std::vec>("double_vector")
            .def(vector_indexing_suite<std::vec>())
            ;

        class_<std::mat>("double_mat")
            .def(vector_indexing_suite<std::mat>())
            ;

        class_<Quad>("Quad", init<int,int>())
            .def("integrate", &Quad::integrate_wrapper)
            .def("integrate_from_string", &Quad::integrate_from_string)
            .def_readonly("nodes", &Quad::nodes)
            .def_readonly("weights", &Quad::weights)
            ;
    }
}
