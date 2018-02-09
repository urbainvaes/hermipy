#include <cmath>
#include <functional>
#include <numeric>
#include <iostream>
#include <tuple>
#include <utility>
#include <vector>

#include <boost/function.hpp>
#include <boost/python.hpp>
#include <boost/math/special_functions/binomial.hpp>

#include "hermite/helpers/templates.hpp"
#include "hermite/helpers/combinatorics.hpp"
#include "hermite/integrator.hpp"
#include "hermite/quadrature.hpp"

#define PI 3.141592653589793238462643383279502884

typedef std::vector< std::vector< std::vector<double> > > cube;
typedef std::vector< std::vector<double> > mat;
typedef std::vector<double> vec;

using namespace std;

namespace hermite {

    // nodes and weight shuold have the right length
    void get_gh_quadrature(int nNodes, vector<double>& nodes, vector<double>& weights) {

        nodes = vec(nNodes);
        weights = vec(nNodes);

        // Unscaled nodes and weights
        vector<double> x;
        vector<double> w;

        switch (nNodes) {
            case 1: nodes = {0.}; weights = {1.}; return; break;
            case 2: x = nodes_2; w = weights_2; break;
            case 4: x = nodes_4; w = weights_4; break;
            case 6: x = nodes_6; w = weights_6; break;
            case 8: x = nodes_8; w = weights_8; break;
            case 10: x = nodes_10; w = weights_10; break;
            case 16: x = nodes_16; w = weights_16; break;
            case 20: x = nodes_20; w = weights_20; break;
            case 30: x = nodes_30; w = weights_30; break;
            case 32: x = nodes_32; w = weights_32; break;
            case 64: x = nodes_64; w = weights_64; break;
            case 100: x = nodes_100; w = weights_100; break;
            default: cout << "Invalid number of nodes (" << nNodes << ") for Gauss-hermite integration" << endl; exit(0);
        }

        // Add symmetric
        for (int i = 0; i < nNodes/2; ++i) {
            nodes[i] = x[i];
            weights[i] = w[i];
            nodes[i + nNodes/2] = -x[i];
            weights[i + nNodes/2] = w[i];
        }

        // Scaling to get rid of factors
        for (int i = 0; i < nNodes; ++i) {
            weights[i] /= sqrt(PI);
            nodes[i] *= sqrt(2);
        }
    }

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

    void get_smolyak_nodes_and_weights(int nVars, vector< vector<double> >& snodes, vector<double>& sweights) {

        snodes = vector< vector<double> > (0);
        sweights = vector<double>(0);

        vector<int> indices = {0, 1, 2, 3, 4, 5, 6};
        vector<int> sizes = {1, 2, 4, 8, 16, 32, 64};

        // Using notations of
        // https://dbarajassolano.wordpress.com/2012/01/26/on-sparse-grid-quadratures/
        int q = indices.size() - 1;

        // Enumeration of multi_indinces
        vector< vector<int> > enum_indices = interval_multi_indices(nVars, q - nVars + 1, q);

        // Auxiliary vectors
        vector< vector<double> > aux_nodes;
        vector<double> aux_weights;

        for (unsigned int i = 0; i < enum_indices.size(); ++i) {

            vector<int> this_index = enum_indices[i];
            vector<int> this_sizes(this_index.size());

            for (unsigned int j = 0; j < this_index.size(); ++j)
                this_sizes[j] = sizes[this_index[j]];

            int sum_index = accumulate(this_index.begin(), this_index.end(), 0.);
            double weight_mod = pow(-1, q - sum_index) * boost::math::binomial_coefficient<double> (nVars - 1, q - sum_index);

            quad_prod(this_sizes, aux_nodes, aux_weights);
            aux_weights = aux_weights * weight_mod;

            snodes.insert(snodes.end(), aux_nodes.begin(), aux_nodes.end());
            sweights.insert(sweights.end(), aux_weights.begin(), aux_weights.end());
        }
    }

    class Quad {
        public:
            Quad(int, int);
            double integrate(boost::function<double(vec)>);
            double integrate_wrapper(boost::python::object);
            mat nodes;
            vec weights;
    };

    Quad::Quad(int nNodes, int nVars) {

        vector<double> nodes_1d(nNodes);
        vector<double> weights_1d(nNodes);

        if (nNodes == 0)
            get_smolyak_nodes_and_weights(nVars, nodes, weights);
        else {
            vector<int> seq(nVars, nNodes);
            quad_prod(seq, nodes, weights);
        }
    }

    double Quad::integrate(boost::function<double(vec)> func) {

        double result = 0.;
        for (unsigned int i = 0; i < nodes.size(); ++i) {
            result += func(nodes[i]) * weights[i];
        }
        return result;
    }

    // ---- PYTHON WRAPPERS ----
    double Quad::integrate_wrapper(boost::python::object func) {
        std::function<double(vec)> lambda;
        switch (this->nodes[0].size()) {
            case 1: lambda = [func](vec v) -> double { return boost::python::extract<double>(func (v[0])); }; break;
            case 2: lambda = [func](vec v) -> double { return boost::python::extract<double>(func(v[0], v[1])); }; break;
            case 3: lambda = [func](vec v) -> double { return boost::python::extract<double>(func(v[0], v[1], v[2])); }; break;
            default: cout << "Dimension must be 1, 2, or 3" << endl; exit(0);
        }
        return integrate(boost::function<double(vec)>(lambda));
    }

    // ---- EXPOSE TO PYTHON ----
    BOOST_PYTHON_MODULE(hermite)
    {
        using namespace boost::python;
        // def("integrate", integrate_wrap);
        // def("integrate_2d", integrate_2d_wrap);
        // def("integrate_3d", integrate_3d_wrap);
        class_<Quad>("Quad", init<int,int>())
            .def("integrate", &Quad::integrate_wrapper)
            ;
    }
}
