#include <iostream>
#include <dlfcn.h>
#include <math.h>

#include "hermite/integrator.hpp"

using namespace std;

double func(vector<double> const& v) {
    // return cos(2*v[0]+2*v[1]+2*v[2]) + v[0]*v[1] + exp(-v[2]*v[2]) + cos(2*v[0]+2*v[1]+2*v[2]) + v[0]*v[1] + exp(-v[2]*v[2]) + cos(2*v[0]+2*v[1]+2*v[2]) + v[0]*v[1] + exp(-v[2]*v[2]) + cos(2*v[0]+2*v[1]+2*v[2]) + v[0]*v[1] + exp(-v[2]*v[2]) + cos(2*v[0]+2*v[1]+2*v[2]) + v[0]*v[1] + exp(-v[2]*v[2]);
    return v[0]*v[0];

}

int main() {
    hermite::Quad quadrature(100, 3);
    double result = quadrature.integrate_from_string("cos(2*v[0]+2*v[1]+2*v[2]) + v[0]*v[1] + exp(-v[2]*v[2]) + cos(2*v[0]+2*v[1]+2*v[2]) + v[0]*v[1] + exp(-v[2]*v[2]) + cos(2*v[0]+2*v[1]+2*v[2]) + v[0]*v[1] + exp(-v[2]*v[2]) + cos(2*v[0]+2*v[1]+2*v[2]) + v[0]*v[1] + exp(-v[2]*v[2]) + cos(2*v[0]+2*v[1]+2*v[2]) + v[0]*v[1] + exp(-v[2]*v[2])");

    cout << "Result = " << result << endl;
}
