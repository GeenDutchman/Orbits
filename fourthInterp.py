#!/usr/bin/python3

_ORDER = 4

def interpolate(x_coords, y_coords):
    if len(y_coords) != len(x_coords):
        raise ValueError('The lengths of x_coords and y_coords must be the same!!')
    if len(y_coords) != _ORDER:
        raise ValueError('Need', _ORDER, 'data points!!')
    def f_dot(x):
        def _skip_term(x, skipped_index):
            product = 1
            for iter in range(_ORDER):
                if iter != skipped_index:
                    product *= x - x_coords[iter]
            return product
        sum = 0
        for i in range(_ORDER):
            sum += (y_coords[i]*_skip_term(x, i))/_skip_term(x_coords[i], i)
        return sum
    return f_dot