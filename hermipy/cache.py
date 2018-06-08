# Copyright (C) 2018 Urbain Vaes
#
# This file is part of hermipy, a python/C++ library for automating the
# Hermite Galerkin method.
#
# hermipy is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# hermipy is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.

import hashlib
import numpy.linalg as la
import sympy as sym
import numpy as np
from functools import wraps
import scipy.sparse as sparse
from hermipy.settings import settings

import os


def gen_hash(extend=None):

    def default_extend(argument):
        raise ValueError("Argument type not supported: " +
                         str(type(argument)))

    if extend is None:
        extend = default_extend

    def the_hash(argument):
        if isinstance(argument, str):
            encoded_str = argument.encode('utf-8')
            return hashlib.md5(encoded_str).hexdigest()
        if isinstance(argument, sparse.csr_matrix):
            a = argument
            new_args = [a.data, a.indices, a.indptr, a.shape]
            return the_hash(new_args)
        if isinstance(argument, np.ndarray):
            return hashlib.md5(argument).hexdigest()
        if isinstance(argument, tuple(sym.core.all_classes)):
            return the_hash(str(argument))
        if isinstance(argument, (list, tuple)):
            hashes = [the_hash(e) for e in argument]
            return the_hash(hash(frozenset(hashes)))
        if isinstance(argument, dict):
            hashes = {kw: the_hash(argument[kw]) for kw in argument}
            return the_hash(hash(frozenset(argument)))
        if isinstance(argument, (int, float)):
            return the_hash(str(hash(argument)))
        if argument is None:
            return the_hash("")
        return the_hash(extend(argument))

    return the_hash


def gen_error(extend=None):

    def default_extend(u, v):
        raise ValueError("Invalid types")

    if extend is None:
        extend = default_extend

    def the_error(u, v):
        if isinstance(u, (float, int)):
            return abs(u - v)
        elif isinstance(u, np.ndarray):
            return la.norm(u - v, 2)
        elif isinstance(u, sparse.csr_matrix):
            assert isinstance(v, sparse.csr_matrix)
            return la.norm(u.data - v.data, 2)
        elif isinstance(u, tuple):
            return sum([the_error(ui, vi) for ui, vi in zip(u, v)])
        return extend(u, v)

    return the_error


def cache(hash_extend=None, error_extend=None, quiet=False):

    hash_fun = gen_hash(extend=hash_extend)
    error_fun = gen_error(extend=error_extend)

    def numpy_extended_load(filename):
        result = np.load(filename)
        if result.shape is ():
            result = float(result)
        return result

    def cache_aux(function):

        @wraps(function)
        def wrapper(*args, **kwargs):
            use_cache = settings['cache']
            if 'cache' in kwargs:
                use_cache = kwargs['cache']
                del kwargs['cache']
            hashes, prefix = [], function.__name__
            for arg in args:
                hashes.append(hash_fun(arg))
            for kw in kwargs:
                hashes.append(hash_fun(kw))
                hashes.append(hash_fun(kwargs[kw]))
            hash_args = hash_fun(('-'.join(hashes)))

            is_sparse = 'sparse' in kwargs and kwargs['sparse'] is True
            ext = '.npz' if is_sparse else '.npy'
            save = sparse.save_npz if is_sparse else np.save
            load = sparse.load_npz if is_sparse else numpy_extended_load

            cachedir = settings['cachedir']
            if not os.path.exists(cachedir):
                os.makedirs(cachedir)

            savefile = cachedir + '/' + prefix + '-' + str(hash_args) + ext
            try:
                result_cache = load(savefile)

            except IOError:
                result = function(*args, **kwargs)
                save(savefile, result)
                return result

            if use_cache:
                return result_cache
            else:
                result = function(*args, **kwargs)
                error = error_fun(result, result_cache)
                if not quiet:
                    assert error < 1e-10
                return result
        return wrapper
    return cache_aux
