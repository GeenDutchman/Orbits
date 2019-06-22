#!/usr/bin/python3

import sys
import numpy as np
from scipy.interpolate import interp1d

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
        second_set[i][0] = second_set[i][0] - 2 * np.pi * cycles - epsilon
    return second_set

def find_residual(extracted, second_set):
    second_set = second_set[second_set[:,0].argsort()]
    f = interp1d(second_set[:,0], second_set[:,1], kind='cubic')
    print('hihihihihihhihihihihihihihi')
    residual = 0
    for i in range(len(extracted)):
        residual += abs(extracted[i][1] - f(second_set[i][0]))
    return residual

def main(epsilon=0.0):
    data = np.loadtxt('binary1.dat', dtype=np.float64)
    
    extracted = orbital_extractor(data)
    if len(extracted) == 0:
        print("Not enough data points!!")
        return 1
    second_set = snipper(extracted, epsilon=float(epsilon))
    residual = find_residual(extracted, second_set)
    print(residual)

if __name__ == "__main__":
    #main(sys.argv[1:])
    if len(sys.argv) > 0:
        main(sys.argv[1])
    else:
        main()
