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

#include <cmath>
#include <dlfcn.h>
#include <fstream>

#include "hermite/types.hpp"
#include "hermite/iterators.hpp"
#include "hermite/discretize.hpp"

using namespace std;

namespace hermite {

void map_point(vec const & translation,
        mat const & dilation,
        double *node,
        double *mapped_node)
{
    u_int dim = translation.size();

    for (u_int i = 0; i < dim; i++)
    {
        mapped_node[i] = translation[i];
        for (u_int j = 0; j < dim; j++)
        {
            mapped_node[i] += dilation[i][j]*node[j];
        }
    }
}

vec discretize(
        s_func func,
        mat const & nodes,
        vec const & translation,
        mat const & dilation)
{
    u_int dim = nodes.size();

    u_int i,j;
    double* node = (double*) malloc(sizeof(double)*dim);
    double* mapped_node = (double*) malloc(sizeof(double)*dim);

    u_int n_points_tot = 1;
    ivec n_points(dim);
    for (u_int i = 0; i < dim; i++)
    {
        n_points[i] = nodes[i].size();
        n_points_tot *= n_points[i];
    }
    vec result(n_points_tot, 0);

    Grid_iterator p(n_points);
    for (i = 0; i < n_points_tot; i++, p.increment())
    {
        for (j = 0; j < dim; j++)
        {
            node[j] = nodes[j][p[j]];
        }
        map_point(translation, dilation, node, mapped_node);
        result[i] = func(mapped_node);
    }

    free(node);
    free(mapped_node);

    return result;
}

void intern_function(string const & function_body)
{
    string name = to_string(hash<string>()(function_body));
    string cpp_file = "/tmp/" + name + ".cpp";
    string so_file = "/tmp/" + name + ".so";
    ifstream test_exists(so_file.c_str());

    if(! test_exists.good()) {
         ofstream write_function;
         write_function.open(cpp_file);
         write_function << "#include <vector>\n#include <cmath>\n";
         write_function << "extern \"C\" double toIntegrate(double *v) {\n";
         write_function << "    return " << function_body << ";\n}";
         write_function.close();

        // Compile file
        string command = "c++ " + cpp_file + " -o " + so_file + " -O3 -Ofast -shared -fPIC";
        system(command.c_str());
    }
}

vec discretize_from_string(
        string function_body,
        mat const & nodes,
        vec const & translation,
        mat const & dilation)
{
    intern_function(function_body);
    string name = to_string(hash<string>()(function_body));
    string so_file = "/tmp/" + name + ".so";
    void *function_so = dlopen(so_file.c_str(), RTLD_NOW);
    s_func func = (s_func) dlsym(function_so, "toIntegrate");
    vec result = discretize(func, nodes, translation, dilation);
    dlclose(function_so);
    return result;
}

}
