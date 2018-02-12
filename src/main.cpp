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

    // typedef double (*integrator_func)(std::vector<double>);
    // typedef void (*voidfunc)();
    // void *myso = dlopen("/home/urbain/scratch/integrator/hello.so", RTLD_NOW);
    // integrator_func func_to_integrate = (integrator_func) dlsym(myso, "hello");
    // voidfunc say_hello = (voidfunc) dlsym(myso, "say_hello");
    // say_hello();
    // result = quadrature.integrate(func_to_integrate);
    // dlclose(myso);

    cout << "Result = " << result << endl;
}
