/*
 * Copyright (C) 2018 Urbain Vaes
 *
 * This file is part of hermipy, a python/C++ library for automating the
 * Hermite Galerkin method.
 *
 * hermipy is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * hermipy is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <https://www.gnu.org/licenses/>.
 */

#include "hermite/templates.hpp"

using namespace std;

hermite::mat operator*(const hermite::mat & mat1, const hermite::mat & mat2) {
    hermite::mat result (mat1.size(), hermite::vec  (mat1.size(), 0.));
    for (unsigned int i = 0; i < mat1.size(); i++) {
        for (unsigned int j = 0; j < mat1.size(); j++) {
            for (unsigned int k = 0; k < mat1.size(); ++k) {
                result[i][j] += mat1[i][k] * mat2[k][j];
            }
        }
    }
    return result;
}

hermite::vec operator*(const hermite::mat & matrix, const hermite::vec & vector) {
    hermite::vec result(matrix.size(),0.);
    for (unsigned int i = 0; i < matrix.size(); i++) {
        for (unsigned int j = 0; j < matrix.size(); j++) {
            result[i] += matrix[i][j] * vector[j];
        }
    }
    return result;
}

double operator*(const hermite::vec & vec1, const hermite::vec & vec2) {
    double result = 0.;
    for (unsigned int i = 0; i < vec1.size(); ++i) {
        result += vec1[i] * vec2[i];
    }
    return result;
}
