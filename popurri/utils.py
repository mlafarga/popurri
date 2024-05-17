"""
General utils
"""
import functools
import time

def time_print(func):
    """
    Decorator to print the time a function takes to run.
    """
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        start_time = time.time()
        result = func(*args, **kwargs)
        total_time = time.time() - start_time
        print(f'{func.__name__} took {total_time} seconds')
        return result
    return wrapper





