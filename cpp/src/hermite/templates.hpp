#ifndef TEMPLATES_H
#define TEMPLATES_H

#include <vector>
#include <cmath>

#include "hermite/types.hpp"


template<class type> std::vector<type> operator-(const std::vector<type>& v1, const std::vector<type>& v2) {
    std::vector<type> result(v1.size());
    for (unsigned int i = 0; i < v1.size(); ++i) {
        result[i] = v1[i] - v2[i];
    }
    return result;
}

template<class type> std::vector<type> operator+(const std::vector<type>& v1, const std::vector<type>& v2) {
    std::vector<type> result(v1.size());
    for (unsigned int i = 0; i < v1.size(); ++i) {
        result[i] = v1[i] + v2[i];
    }
    return result;
}

template<class type> double fabs(std::vector<type> vec) {
    double result = 0.;
    for (unsigned int i = 0; i < vec.size(); ++i) {
        double r = fabs(vec[i]);
        result += r*r;
    }
    return sqrt(result);
}

template<class type> std::vector<type> operator*(const std::vector<type>& vec, const double& x) {
    std::vector<type> result = vec;
    for (unsigned int i = 0; i < vec.size(); ++i) {
        result[i] = result[i]*x;
    }
    return result;
}

hermite::mat operator*(const hermite::mat &, const hermite::mat &);
hermite::vec operator*(const hermite::mat &, const hermite::vec &);
double operator*(const hermite::vec &, const hermite::vec &);

#endif
