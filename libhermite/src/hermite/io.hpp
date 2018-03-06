#ifndef IO_H
#define IO_H

#include "types.hpp"
#include <string>

namespace hermite {

void printVec(std::vec);
void printMat(std::mat, std::string s = "  ");
void printCube(std::cube, std::string s = " ");

}

#endif
