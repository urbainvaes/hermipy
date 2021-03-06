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

from functools import wraps
import hermipy.settings as rc
import numpy as np
import time

stats = {}
indent = 0
threshold = 0


def debug(force=False):
    def decorator(function):
        def stringify_arg(arg):
            if isinstance(arg, list):
                return str([stringify_arg(a) for a in arg])
            if isinstance(arg, np.ndarray):
                return "numpy ndarray " + str(arg.shape)
            else:
                return str(arg)

        @wraps(function)
        def wrapper(*args, **kwargs):
            key = function.__name__ + '-' + function.__module__
            if rc.settings['debug'] or force:
                print("▶ Entering function " + key)
                print("-▲ args")
                for a in args:
                    print("--● A " + stringify_arg(a))
                print("-▲ kwargs")
                for kw in kwargs:
                    print("--● KW " + str(kw) + ": " +
                          stringify_arg(kwargs[kw]))
            result = function(*args, **kwargs)
            return result
        return wrapper
    return decorator


def log_stats(force=False):
    def decorator(function):
        @wraps(function)
        def wrapper(*args, **kwargs):
            global indent, threshold

            indent_str = ''.join([" "]*indent) + str(indent) + " "
            key = function.__name__ + '-' + function.__module__
            if rc.settings['trails'] or force:
                print(indent_str + "Entering " + key)

            indent += 1
            time_start = time.time()
            result = function(*args, **kwargs)
            time_end = time.time()
            indent -= 1

            if key not in stats:
                stats[key] = {'Calls': 0, 'Time': 0}
            stats[key]['Calls'] += 1
            stats[key]['Time'] += time_end - time_start
            if rc.settings['trails'] or force:
                delta_time = time_end - time_start
                if delta_time > threshold:
                    print(indent_str + "Function " + key
                          + ": " + str(delta_time))
            return result
        return wrapper
    return decorator


def print_stats():
    print('\n-- [statistics] --')
    for key in stats:
        _calls, _time = str(stats[key]['Calls']), str(stats[key]['Time'])
        print("| '" + key + "': Calls: " + _calls + ", Time: " + _time)
    print('--')
