#include <unordered_map>
#include <boost/math/special_functions/binomial.hpp>
#include "hermite/iterators.hpp"
#include "hermite/lib.hpp"
#include "hermite/types.hpp"

using namespace std;
using boost::math::binomial_coefficient;

namespace hermite
{

u_int bissect_degree(u_int dim, u_int n_polys)
{
    u_int degree_1 = 0, degree_2 = 150;

    int img_1 = (int) binomial_coefficient<double> (degree_1 + dim, dim) - (int) n_polys;
    int img_2 = (int) binomial_coefficient<double> (degree_2 + dim, dim) - (int) n_polys;

    if (img_1 > 0 || img_2 < 0)
    {
        cout << "Can't find degree, Invalid arguments!" << endl;
        exit(0);
    }

    if (img_1 == 0)
        return degree_1;

    if (img_2 == 0)
        return degree_2;

    while (true)
    {
        u_int new_degree = (degree_1 + degree_2)/2;
        int new_img = (int) binomial_coefficient<double> (new_degree + dim, dim) - (int) n_polys;

        if (new_img < 0)
        {
            degree_1 = new_degree;
            img_1 = new_img;
        }
        else if (new_img > 0)
        {
            degree_2 = new_degree;
            img_2 = new_img;
        }
        else
        {
            return new_degree;
        }
    }
}

u_int hash_multi_ind(ivec v, int degree)
{
    u_int base = degree + 1;
    u_int result = 0;
    u_int unit = 1;
    for(u_int i = 0; i < v.size(); i++)
    {
        result += v[i]*unit;
        unit *= base;
    }
    return result;
}

unordered_map<u_int, u_int> hash_table(u_int dim, u_int degree)
{
    Multi_index_iterator m(dim, degree);
    unordered_map<u_int, u_int> lin_indices;
    for (u_int i = 0; !m.isFull(); i++, m.increment())
    {
        u_int hash = hash_multi_ind(m.get(), degree);
        lin_indices.insert(pair<u_int, u_int>(hash, i));
    }
    return lin_indices;
}

bool isAligned(const ivec & m, u_int dir)
{
    for (u_int j = 0; j < m.size(); j++)
    {
        if (m[j] != 0 && j != dir)
        {
            return false;
        }
    }
    return true;
}

bool isAligned(const ivec & m, const ivec & dirs)
{
    for (u_int j = 0; j < m.size(); j++)
    {
        bool in_dirs = false;
        for(u_int k = 0; k < dirs.size(); k++)
        {
            if (j == dirs[k])
            {
                in_dirs = true;
            }
        }

        if (m[j] != 0 && !in_dirs)
        {
            return false;
        }
    }
    return true;
}

}
