/*! \file integrator.hpp
 * Functions and classes related to quadratures
 */

#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include <string>

#include <boost/function.hpp>
#include <boost/python.hpp>

namespace hermite {
    class Quad {
        public:
            Quad(int, int);
            double integrate(boost::function<double(std::vector<double> const&)> const&);
            double integrate_from_string(std::string const& function_body);
            double integrate_wrapper(boost::python::object const&);
            std::vector< std::vector<double> > nodes;
            std::vector<double> weights;
    };
}

#endif
