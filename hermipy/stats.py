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

import hermipy.settings as rc
from functools import wraps
import time

stats = {}
indent = 0
threshold = 1e-3


def log_stats(function):
    @wraps(function)
    def wrapper(*args, **kwargs):
        global indent, threshold
        indent += 1
        time_start = time.time()
        result = function(*args, **kwargs)
        time_end = time.time()
        indent -= 1

        key = function.__name__ + '-' + function.__module__
        if key not in stats:
            stats[key] = {'Calls': 0, 'Time': 0}
        stats[key]['Calls'] += 1
        stats[key]['Time'] += time_end - time_start
        if rc.settings['trails']:
            delta_time = time_end - time_start
            if delta_time > threshold:
                print(str(indent) + " Function " + key
                      + ": " + str(delta_time))
        return result
    return wrapper


def print_stats():
    print('\n-- [statistics] --')
    for key in stats:
        _calls, _time = str(stats[key]['Calls']), str(stats[key]['Time'])
        print("| '" + key + "': Calls: " + _calls + ", Time: " + _time)
    print('--')
