from functools import wraps
from .settings import settings
import time

stats = {}


def log_stats(function):
    @wraps(function)
    def wrapper(*args, **kwargs):
        time_start = time.time()
        result = function(*args, **kwargs)
        time_end = time.time()

        key = function.__name__
        if key not in stats:
            stats[key] = {'Calls': 0, 'Time': 0}
        stats[key]['Calls'] += 1
        stats[key]['Time'] += time_end - time_start
        if settings['trails']:
            print("Function " + key + ": " + str(time_end - time_start))
        return result
    return wrapper
