#ifndef HERMITE_ITERATORS_H
#define HERMITE_ITERATORS_H

#include <string>
#include "hermite/types.hpp"

namespace hermite {

class Vector_iterator {

    protected:

    unsigned int dim;
    std::ivec multi_index;
    bool full;

    public:

    const std::ivec& get() const {
        return multi_index;
    }

    int operator[](int i) {
        return multi_index[i];
    }

    std::ivec get() {
        return multi_index;
    }

    bool isFull() const {
        return full;
    }

    public:

    virtual void increment() = 0;
    virtual void reset() = 0;
    Vector_iterator(int dim): dim(dim), multi_index(std::ivec(dim, 0)), full(false) {}

};

class Multi_index_iterator : public Vector_iterator {

    unsigned int sum;

    // Upper bound included (like polynomial degree)
    const unsigned int upper_bound;

    public:

    void increment();
    void reset();
    Multi_index_iterator(unsigned int dim, unsigned int upper_bound):
        Vector_iterator(dim), sum(0), upper_bound(upper_bound) {}
};

class Hyper_cube_iterator : public Vector_iterator {

    // Upper bounds excluded
    const std::ivec upper_bounds;

    public:

    void increment();
    void reset();
    Hyper_cube_iterator(const std::ivec & upper_bounds):
        Vector_iterator(upper_bounds.size()), upper_bounds(upper_bounds) {}
};

std::imat list_multi_indices(u_int dim, u_int upper_bound);
std::imat list_cube_indices(const std::ivec & upper_bounds);

}

#endif
