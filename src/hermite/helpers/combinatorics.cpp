#include "hermite/helpers/combinatorics.hpp"

#include <iostream>
#include <math.h>

using namespace std;

// Standard normal gaussian
double gaussian(vector<double> y) {
    double result = 1.;
    for (unsigned int i = 0; i < y.size(); ++i) {
        result *= exp(-y[i]*y[i]/2)/(sqrt(2*M_PI));
    }
    return result;
}

vector< vector<int> > equal_multi_indices(int n, int c) {

    // Univariate case
    if (n == 1)
        return vector< vector<int> > (1, vector<int> (1, c));

    // Case where |i| = 0
    if (c == 0)
        return vector< vector<int> > (1, vector<int> (n, 0));

    // Initialization of result
    vector< vector<int> > result(0);

    // Iteration on the first variable
    for (int i = 0; i <= c; ++i) {
        vector< vector<int> > aux = equal_multi_indices(n-1, c-i);
        for (unsigned int j = 0; j < aux.size(); ++j) {
            aux[j].push_back(i);
            result.push_back(aux[j]);
        }
    }

    return result;
}

vector< vector<int> > interval_multi_indices(int n, int a, int b) {

    if (b == a)
        return equal_multi_indices(n, a);

    // Auxiliary results
    vector< vector<int> > r_lower = interval_multi_indices(n, a, b-1);
    vector< vector<int> > r_equal = equal_multi_indices(n, b);

    // Concatenation of the vectors
    r_lower.insert(r_lower.end(), r_equal.begin(), r_equal.end());

    return r_lower;
}

vector< vector<int> > lower_multi_indices(int n, int b) {
    return interval_multi_indices(n, 0, b);
}
