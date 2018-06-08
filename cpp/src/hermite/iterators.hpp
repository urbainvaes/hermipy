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

#ifndef HERMITE_ITERATORS_H
#define HERMITE_ITERATORS_H

#include <string>
#include "hermite/types.hpp"

namespace hermite {

class Vector_iterator {

    protected:

        unsigned int const dim;
        ivec multi_index;
        bool full;

    public:

        const ivec & get() const {
            return multi_index;
        }

        int operator[](int i) const {
            return multi_index[i];
        }

        bool isFull() const {
            return full;
        }

    public:

        virtual void increment() = 0;
        virtual void reset() = 0;
        Vector_iterator(int dim): dim(dim), multi_index(ivec(dim, 0)), full(false) {}

};

class Multi_index_iterator : public Vector_iterator {

    unsigned int sum;

    // Upper bound included (like polynomial degree)
    const unsigned int upper_bound;

    public:

    // Get linear index from multi-index
    static u_int index(const ivec & m_vec);
    static u_int size(u_int degree, u_int dim);

    unsigned int get_sum()
    {
        return sum;
    }

    void increment();
    void reset();
    Multi_index_iterator(unsigned int dim, unsigned int upper_bound):
        Vector_iterator(dim), sum(0), upper_bound(upper_bound) {}
};

class Hyper_cube_iterator : public Vector_iterator {

    // Upper bounds excluded
    const ivec upper_bounds;

    public:

    void increment();
    void reset();
    Hyper_cube_iterator(const ivec & upper_bounds):
        Vector_iterator(upper_bounds.size()), upper_bounds(upper_bounds) {}
};

imat list_multi_indices(u_int dim, u_int upper_bound);
imat list_cube_indices(const ivec & upper_bounds);

}

#endif
