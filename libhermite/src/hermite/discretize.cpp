#include <cmath>

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

    Hyper_cube_iterator p(n_points);
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

}
