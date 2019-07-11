#!/usr/bin/python3

import numpy as np
from find_crossing import crossing_finder
from fourthInterp import interpolate, _ORDER
from scipy.optimize import minimize

data = np.genfromtxt('binary1.dat', dtype=np.float64, names=True)
time_col = data['time']
angle_col = data['star_angle']
r_dot_col = data['star_r_dot']

indecies, t_crossings, r_crossings = crossing_finder(time_col, r_dot_col)

prev_angle = 0
for i in range(len(t_crossings)):
    time_to_pos = interpolate(t_crossings[i], r_crossings[i])
    def abs_of_ttp(x):
        return abs(time_to_pos(x))
    result = minimize(abs_of_ttp, np.average(t_crossings[i][1:3]))

    cross_index = indecies[i]
    time_to_angle = interpolate(t_crossings[i], angle_col[cross_index:cross_index + _ORDER])
    cumulative_angle = time_to_angle(result.x[0])
    #print('At time:', result.x[0], 'Accumulative angle:', cumulative_angle, 'Single Orbit:', cumulative_angle - prev_angle)
    print(t_crossings[i], r_crossings[i],
          angle_col[cross_index:cross_index + _ORDER])
    prev_angle = cumulative_angle


