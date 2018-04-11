#include <cmath>

#include "hermite/hermite.hpp"

using namespace std;

namespace hermite {

void hermite_eval(double x,
        u_int degree,
        vec & values)
{

    values[0] = 1;
    values[1] = x;

    for (unsigned i = 1; i < degree; i++)
    {
        values[i+1] = REC_A(i)*x*values[i] - REC_B(i)*values[i-1];
    }
}}
