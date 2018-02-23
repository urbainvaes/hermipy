#include "hermite/templates.hpp"

using namespace std;

mat operator*(const mat& mat1, const mat& mat2) {
    mat result (mat1.size(), vec (mat1.size(), 0.));
    for (unsigned int i = 0; i < mat1.size(); i++) {
        for (unsigned int j = 0; j < mat1.size(); j++) {
            for (unsigned int k = 0; k < mat1.size(); ++k) {
                result[i][j] += mat1[i][k] * mat2[k][j];
            }
        }
    }
    return result;
}

vec operator*(const mat& matrix, const vec& vector) {
    vec result(matrix.size(),0.);
    for (unsigned int i = 0; i < matrix.size(); i++) {
        for (unsigned int j = 0; j < matrix.size(); j++) {
            result[i] += matrix[i][j] * vector[j];
        }
    }
    return result;
}

double operator*(const vec& vec1, const vec& vec2) {
    double result = 0.;
    for (unsigned int i = 0; i < vec1.size(); ++i) {
        result += vec1[i] * vec2[i];
    }
    return result;
}
