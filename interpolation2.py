#!/usr/bin/python3

import numpy as np
from scipy.interpolate import interp1d

data = np.loadtxt('binary1.dat', dtype=np.float64)

def orbital_extractor(data, cycles=1):
    extractor = []
    for i in range(len(data)):
        if data[i][-1] > 3 * np.pi * cycles:
            break
        elif data[i][-1] >= np.pi * cycles:
            extractor.append([data[i][-1], data[i][-2]])
            extracted = np.array(extractor)
    return extracted

def snipper(extracted, cycles=1, epsilon=0.0):
    second_set = np.copy(extracted)
    for i in range(len(second_set)):
        second_set[i][0] -= 2 * np.pi * cycles - epsilon
    return second_set

def find_residual(extracted, second_set):
    f = interp1d(second_set[:,0], second_set[:,1], kind='cubic')
    residual = 0
    for i in range(len(extracted)):
        residual += abs(extracted[i][1] - f(extracted[i][1]))
    return residual
