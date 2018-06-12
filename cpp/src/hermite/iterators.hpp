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

namespace hermite
{

class Vector_iterator
{
    protected:

        u_int const dim;
        ivec multi_index;
        bool full;

    public:

        const ivec & get() const
        {
            return multi_index;
        }

        int operator[](int i) const
        {
            return multi_index[i];
        }

        bool isFull() const
        {
            return full;
        }

        void reset()
        {
            full = false;
            for (unsigned int i = 0; i < dim; i++)
            {
                multi_index[i] = 0;
            }
        }

        virtual void increment() = 0;
        Vector_iterator(int dim): dim(dim), multi_index(ivec(dim, 0)), full(false) {}

};

class Hyperbolic_cross_iterator : public Vector_iterator
{
    // Upper bound included (like polynomial degree)
    const u_int upper_bound;
    const imat list;
    u_int index;

    public:

    static imat list(u_int dim, u_int upper_bound);
    void increment();
    Hyperbolic_cross_iterator(u_int dim, u_int upper_bound):
        Vector_iterator(dim), upper_bound(upper_bound) {}
};

class Multi_index_iterator : public Vector_iterator
{
    // Upper bound included (like polynomial degree)
    const u_int upper_bound;

    public:

    // Get linear index from multi-index
    static u_int index(const ivec & m_vec);
    static u_int size(u_int degree, u_int dim);
    static imat list(u_int dim, u_int upper_bound);

    void increment();
    Multi_index_iterator(u_int dim, u_int upper_bound):
        Vector_iterator(dim), upper_bound(upper_bound) {}
};

class Hyper_cube_iterator : public Vector_iterator
{
    // Upper bounds excluded
    const ivec upper_bounds;

    public:

    static imat list(const ivec & upper_bounds);

    void increment();
    Hyper_cube_iterator(const ivec & upper_bounds):
        Vector_iterator(upper_bounds.size()), upper_bounds(upper_bounds) {}
};

}

#endif
