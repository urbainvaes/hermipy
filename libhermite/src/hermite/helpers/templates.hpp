#ifndef TEMPLATES_H
#define TEMPLATES_H

#include <vector>
#include <cmath>

namespace std {
    typedef std::vector< std::vector< std::vector<double> > > cube;
    typedef std::vector< std::vector<double> > mat;
    typedef std::vector<double> vec;
}

/*!
 * Compute the difference between two vectors. Using this function requires that the difference between elements of v1 and v2 be well-defined.
 */
template<class type> std::vector<type> operator-(const std::vector<type>& v1, const std::vector<type>& v2) {
    std::vector<type> result(v1.size());
    for (unsigned int i = 0; i < v1.size(); ++i) {
        result[i] = v1[i] - v2[i];
    }
    return result;
}

/*!
 * Compute the sum of two vectors. Using this function requires that the sum of elements of v1 and v2 be well-defined.
 */
template<class type> std::vector<type> operator+(const std::vector<type>& v1, const std::vector<type>& v2) {
    std::vector<type> result(v1.size());
    for (unsigned int i = 0; i < v1.size(); ++i) {
        result[i] = v1[i] + v2[i];
    }
    return result;
}

/*!
 * Compute the absolute value of a vector, by computing the absolute value of individual elements.
 */ 
template<class type> double fabs(std::vector<type> vec) {
    double result = 0.;
    for (unsigned int i = 0; i < vec.size(); ++i) {
        double r = fabs(vec[i]);
        result += r*r;
    }
    return sqrt(result);
}

/*!
 * Compute the product between a vector and a double.
 */ 
template<class type> std::vector<type> operator*(const std::vector<type>& vec, const double& x) {
    std::vector<type> result = vec;
    for (unsigned int i = 0; i < vec.size(); ++i) {
        result[i] = result[i]*x;
    }
    return result;
}

/*!
 * Multiply matrix 1 by matrix 2
 */
std::mat operator*(const std::mat& mat1, const std::mat& mat2);
std::vec operator*(const std::mat& mat, const std::vec& vec);
double operator*(const std::vec& vec1, const std::vec& vec2);

#endif
