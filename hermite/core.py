import pdb

import scipy.sparse as ss
import numpy as np
import hermite_cpp as hm

from hermite.settings import settings
from hermite.stats import log_stats
from hermite.cache import cache


# Conversion functions {{{


def convert_to_cpp_vec(vec):
    cpp_vec = hm.double_vec()
    cpp_vec.extend(vec)
    return cpp_vec


def convert_to_cpp_mat(mat):
    cpp_mat = hm.double_mat()
    for vec in mat:
        cpp_mat.append(convert_to_cpp_vec(vec))
    return cpp_mat


def convert_to_cpp_cube(cube):
    cpp_cube = hm.double_cube()
    for mat in cube:
        cpp_cube.append(convert_to_cpp_mat(mat))
    return cpp_cube


@log_stats
def to_cpp_array(*args):
    if len(args) > 1:
        return (to_cpp_array(arg) for arg in args)
    array, dim = args[0], 0
    if type(array) in (list, np.ndarray):
        dim = 1
        if type(array[0]) in (list, np.ndarray):
            dim = 2
            if type(array[0][0]) in (list, np.ndarray):
                dim = 3
    if dim is 1:
        array = convert_to_cpp_vec(array)
    elif dim is 2:
        array = convert_to_cpp_mat(array)
    elif dim is 3:
        array = convert_to_cpp_cube(array)
    return array


@log_stats
def convert_to_cpp_sparse(mat):
    assert type(mat) is ss.csr_matrix
    data = hm.double_vec()
    indices = hm.int_vec()
    indptr = hm.int_vec()
    data.extend(mat.data)
    indices.extend(mat.indices.astype('uint32'))
    indptr.extend(mat.indptr.astype('uint32'))
    size1, size2 = mat.shape[0], mat.shape[1]
    return hm.to_spmat(data, indices, indptr, size1, size2)


@log_stats
def to_csr(sp_matrix):
    assert isinstance(sp_matrix, hm.sparse_matrix)
    rcv = np.array(hm.row_col_val(sp_matrix))
    shape = (sp_matrix.size1(), sp_matrix.size2())
    return ss.csr_matrix((rcv[2], (rcv[0], rcv[1])), shape=shape)


@log_stats
def to_cpp(arg):
    if type(arg) is np.ndarray or type(arg) is list:
        return to_cpp_array(arg)
    elif type(arg) is ss.csr_matrix:
        return convert_to_cpp_sparse(arg)
    else:
        return arg


@log_stats
def to_numpy(arg):
    if type(arg) == hm.double_vec:
        return np.array(arg)
    if type(arg) == hm.double_mat:
        return hm.to_numpy(arg)
    elif type(arg) == hm.double_cube:
        return np.array(arg)
    elif type(arg) == hm.sparse_matrix:
        return to_csr(arg)


def to_numeric(var):

    if isinstance(var, list):
        return [to_numeric(v) for v in var]

    var_letters = {'x': 0, 'y': 1, 'z': 2}
    var_array = {'v[' + str(i) + ']': i for i in range(10)}
    var_dict = {**var_letters, **var_array}
    if var in var_dict:
        return var_dict[var]
    else:
        return var
# }}}
# Wrapper functions {{{


@cache()
@log_stats
def discretize(function, nodes, translation, dilation):
    nodes, translation, dilation = to_cpp_array(nodes, translation, dilation)
    return np.array(hm.discretize(function, nodes, translation, dilation))


@cache()
@log_stats
def integrate(fgrid, nodes, weights):
    fgrid, nodes, weights = to_cpp_array(fgrid, nodes, weights)
    result = hm.integrate(fgrid, nodes, weights)
    return result


@cache()
@log_stats
def transform(degree, fgrid, nodes, weights, forward):
    fgrid, nodes, weights = to_cpp_array(fgrid, nodes, weights)
    return np.array(hm.transform(degree, fgrid, nodes, weights, forward))


@cache()
@log_stats
def triple_products(degree):
    return np.array(hm.triple_products(degree))


@cache()
@log_stats
def varf(degree, fgrid, nodes, weights, sparse=False):
    if settings['debug']:
        print("Entering varf with dim: " + str(len(nodes)))
    fgrid, nodes, weights = to_cpp_array(fgrid, nodes, weights)
    args = [degree, fgrid, nodes, weights]
    function = hm.varf_sp if sparse else hm.varf
    return log_stats(to_numpy)(log_stats(function)(*args))


@cache()
@log_stats
def varfd(dim, degree, direction, var, numpy=True):
    var = to_cpp(var)
    result = log_stats(hm.varfd)(dim, degree, direction, var)
    return log_stats(to_numpy)(result)


@log_stats
def tensorize(inp, dim=None, direction=None, sparse=False):
    is_scalar = isinstance(inp[0], (float, int))
    if is_scalar and dim is None and direction is None:
        return np.prod(inp)
    if isinstance(inp[0], ss.csr_matrix):
        inp = [np.array(m.todense()) for m in inp]
    inp = to_cpp_array(inp)
    if dim is not None and direction is not None:
        args = [inp, dim, direction]
    else:
        args = [inp]
    if sparse:
        tensorize_fun = hm.tensorize_sp
        convert_fun = to_csr
    else:
        tensorize_fun = hm.tensorize
        convert_fun = np.array
    return convert_fun(tensorize_fun(*args))


@cache()
@log_stats
def project(inp, dim, direction):
    inp = to_cpp_array(inp)
    direction = to_numeric(direction)
    return np.array(hm.project(inp, dim, direction))


@cache()
@log_stats
def multi_indices(dim, degree):
    result = hm.list_multi_indices(dim, degree)
    return np.asarray(result, dtype=int)

