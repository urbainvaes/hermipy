#ifndef INNER_H
#define INNER_H

#include <assert.h>
#include <iostream>
#include "hermite/types.hpp"

namespace hermite
{

vec inner(const vec & s1,
          const vec & s2,
          const ivec & dirs1,
          const ivec & dirs2);
}

#endif
