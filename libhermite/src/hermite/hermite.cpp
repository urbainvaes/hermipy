#include <cmath>

#include "hermite/hermite.hpp"

using namespace std;

namespace hermite {

void hermite_eval(double x,
        u_int degree,
        vec & values,
        bool l2)
{

    double factor = l2 ? sqrt(1 / sqrt(2*M_PI) *  exp(-x*x/2)) : 1;
    // cout << factor << endl;
    values[0] = 1*factor;
    values[1] = x*factor;

    for (unsigned i = 1; i < degree; i++)
    {
        values[i+1] = REC_A(i)*x*values[i] - REC_B(i)*values[i-1];
    }
}}
