from libhermite import hermite as hm
import numpy as np
import inspect


def convert_to_cpp_vec(vec):
    cpp_vec = hm.double_vec()
    cpp_vec.extend(vec)
    return cpp_vec


def convert_to_cpp_mat(mat):
    cpp_mat = hm.double_mat()
    for vec in mat:
        cpp_mat.append(convert_to_cpp_vec(vec))
    return cpp_mat


def convert_to_cpp_array(array):
    dim = 0
    if type(array) in (list, np.ndarray):
        dim = 1
        if type(array[0]) in (list, np.ndarray):
            dim = 2
    if dim is 1:
        array = convert_to_cpp_vec(array)
    elif dim is 2:
        array = convert_to_cpp_mat(array)
    return array


def convert_to_cpp(*names):
    def convert(function):
        sig = inspect.signature(function)

        def wrapper(*args, **kwargs):
            ba = sig.bind(*args, **kwargs)
            args_od = ba.arguments
            for key in args_od:
                if key in names:
                    arg = args_od[key]
                    args_od[key] = convert_to_cpp_array(arg)
            return function(**args_od)
        return wrapper
    return convert


@convert_to_cpp('nodes', 'translation', 'dilation')
def discretize(function, nodes, translation, dilation):
    return np.array(hm.discretize(function, nodes, translation, dilation))


@convert_to_cpp('fgrid', 'nodes', 'weights')
def integrate(fgrid, nodes, weights):
    return hm.integrate(fgrid, nodes, weights)


@convert_to_cpp('fgrid', 'nodes', 'weights')
def transform(degree, fgrid, nodes, weights, forward):
    return np.array(hm.transform(degree, fgrid, nodes, weights, forward))


def triple_products(degree):
    return np.array(hm.triple_products(degree))


@convert_to_cpp('fgrid', 'nodes', 'weights')
def varf(degree, fgrid, nodes, weights):
    return np.array(hm.varf(degree, fgrid, nodes, weights))


@convert_to_cpp('var')
def dvarf(dim, degree, direction, var):
    return np.array(hm.dvarf(dim, degree, direction, var))
