/*! \file combinatorics.hpp
 * Helper functions for combinatorics operations
 */ 

#ifndef COMBINATORIAL_H
#define COMBINATORIAL_H

#include <vector>

/*!
 * Calculate all the n-dimensional multi-indices \f$ \alpha \f$ such that \f$ |\alpha| = c \f$.
 * \return A vector whose elements are the multi-indices.
 */ 
std::vector< std::vector<int> > equal_multi_indices(int n, int c);

/*!
 * Convenience function to get all the n-dimensional multi-indices \f$ \alpha \f$ such that \f$ a \leq |\alpha| \leq b \f$.
 */ 
std::vector< std::vector<int> > interval_multi_indices(int n, int a, int b);

/*!
 * Convenience function to get all the n-dimensional multi-indices \f$ \alpha \f$ such that \f$ |\alpha| \leq b \f$.
 */ 
std::vector< std::vector<int> > lower_multi_indices(int n, int b);


/*!
 * Evaluate the standard gaussian probability density at \f$ y \f$.
 * \f[
 *      \frac{1}{2\pi} \exp(- \sum_{i=1}^n y_i^2 )
 * \f]
 */ 
double gaussian(std::vector<double> y);

#endif
