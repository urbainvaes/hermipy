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
#include "hermite/io.hpp"
#include "hermite/lib.hpp"
#include "hermite/types.hpp"
#include <unordered_map>

namespace hermite
{

class Vector_iterator
{
    protected:
        u_int const dim;
        ivec multi_index;
        bool full;

        virtual ~Vector_iterator() = default;
        Vector_iterator(int dim = 0):
            dim(dim), multi_index(ivec(dim, 0)), full(false) {}

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

        virtual void reset()
        {
            full = false;
            for (unsigned int i = 0; i < dim; i++)
            {
                multi_index[i] = 0;
            }
        }

        virtual void increment() = 0;
};

// Iterator on grids
class Grid_iterator : public Vector_iterator
{
    private:
        // Upper bounds excluded
        const ivec upper_bounds;

    public:
        static imat list(const ivec & upper_bounds);

        // Get linear index from multi-index
        u_int index(const ivec & m_vec);

        void increment();
        Grid_iterator(const ivec & upper_bounds):
            Vector_iterator(upper_bounds.size()), upper_bounds(upper_bounds) {}
};

// Iterators on multi-indices
class Multi_index_iterator : public Vector_iterator
{
    protected:
        std::unordered_map<std::string,u_int> hash_table;
        u_int index_list;
        imat list;

        Multi_index_iterator(u_int dim): Vector_iterator(dim), index_list(0) {}

    public:
        virtual ~Multi_index_iterator() = default;

        imat getList()
        {
            return list;
        }

        u_int size()
        {
            return list.size();
        }

        bool has(const ivec & m_vec)
        {
             std::string hash_mvec = hash_print(m_vec);
             return ! (hash_table.find(hash_mvec) == hash_table.end());
        }

        void reset()
        {
            Vector_iterator::reset();
            index_list = 0;
        }

        u_int index(const ivec & m_vec);
        void increment();
};

// Consistent orderings
class Cube_iterator : public Multi_index_iterator
{
    private:
        const u_int degree;

    public:
        Cube_iterator(u_int dim, u_int degree);

        static bool s_increment(ivec & multi_index, u_int degree);
        static imat s_list(u_int dim, u_int degree);
        static u_int s_index(const ivec & m_vec);
        static u_int s_size(u_int dim, u_int degree);
        static u_int s_get_degree(u_int dim, u_int n_polys);
};

class Triangle_iterator : public Multi_index_iterator
{
    private:
        const u_int degree;

    public:
        Triangle_iterator(u_int dim, u_int degree);

        static bool s_increment(ivec & multi_index, u_int degree);
        static imat s_list(u_int dim, u_int degree);
        static u_int s_index(const ivec & m_vec);
        static u_int s_size(u_int dim, u_int degree);
        static u_int s_get_degree(u_int dim, u_int n_polys);
        static u_int s_find_dim(u_int degree, u_int n_polys);

};

class Cross_iterator : public Multi_index_iterator
{
    private:
        const u_int degree;

    public:
        Cross_iterator(u_int dim, u_int degree);

        static bool s_increment(ivec & multi_index, u_int degree);
        static imat s_list(u_int dim, u_int degree);
        static u_int s_size(u_int dim, u_int degree);
        static u_int s_get_degree(u_int dim, u_int n_polys);
};

// Non consistent (subdegree wouldn't work)
class Cross_iterator_nc : public Multi_index_iterator
{
    private:
        const u_int degree;

    public:
        Cross_iterator_nc(u_int dim, u_int degree);

        static bool s_increment(ivec & multi_index, u_int degree);
        static imat s_list(u_int dim, u_int degree);
        static u_int s_size(u_int dim, u_int degree);
        static u_int s_get_degree(u_int dim, u_int n_polys);
};

class Rectangle_iterator : public Multi_index_iterator
{
    private:
        const u_int degree;

    public:
        Rectangle_iterator(u_int dim, u_int degree);

        static bool s_increment(ivec & multi_index, u_int degree);
        static imat s_list(u_int dim, u_int degree);
        static u_int s_size(u_int dim, u_int degree);
        static u_int s_get_degree(u_int dim, u_int n_polys);
};

}

#endif
