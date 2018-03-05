#ifndef HERMITE_HERMITE_H
#define HERMITE_HERMITE_H

#include "hermite/types.hpp"

// Coefficients in recurrence Hermite polynomials
#define REC_A(i) (1/sqrt(i+1))
#define REC_B(i) (sqrt(i)/sqrt(i+1))

namespace hermite {

void hermite_eval(double x,
        std::u_int degree,
        std::vec & values,
        bool l2);

}

#endif
