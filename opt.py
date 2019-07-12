#!/usr/bin/python3
import sys
import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import minimize


def find_residual(phi_old, r_old, phi_new, r_new):  # new is broader than old
    residuals = []
    func = interp1d(phi_new, r_new, kind='cubic')
    for i in range(len(phi_old)):
        r_interp = func(phi_old[i])
        residuals.append(r_old[i] - r_interp)
    return np.linalg.norm(residuals)


def shifter(phi_column, epsilon, nOrbits):
    phi_shifted = [x - (2*np.pi*nOrbits + epsilon) for x in phi_column]
    return phi_shifted


def copy_piece(phi_column, r_column, start_cut, end_cut):
    phi = []
    r = []
    for row in range(len(phi_column)):
        if phi_column[row] > end_cut:
            break
        elif phi_column[row] >= start_cut:
            phi.append(phi_column[row])
            r.append(r_column[row])
    return phi, r


def opt(epsilon, *args):
    n = args[0]
    epsilon = epsilon[0]
    old_lower_bound = 2 * n * np.pi - np.pi
    old_upper_bound = 2 * n * np.pi + np.pi
    #print(old_lower_bound, old_upper_bound)
    phi_old, r_old = copy_piece(
        phi_original, r_original, old_lower_bound, old_upper_bound)
    new_lower_bound = 2 * n * np.pi - 3/2 * np.pi
    new_upper_bound = 2 * n * np.pi + 3/2 * np.pi
    #print(new_lower_bound, new_upper_bound)
    phi_shifted = shifter(phi_original, epsilon, 1)
    phi_new, r_new = copy_piece(
        phi_shifted, r_original, new_lower_bound, new_upper_bound)
    #print(epsilon)

    # if phi_old[-1] >= phi_new[-1]:
    #     raise ValueError('The data needs to have more orbits!!')
    return find_residual(phi_old, r_old, phi_new, r_new)


i = 1

file_name = "R100.dat"
while i < len(sys.argv):
    if sys.argv[i] == "--read" or sys.argv[i] == "-r":
        i += 1
        file_name = sys.argv[i]
    else:
        print('\n "', sys.argv[i], '" is not an option!!')
        exit(1)
    i += 1    
data = np.genfromtxt(file_name, dtype=np.float64, names=True)
phi_original = data['star_angle']
r_original = data['star_r']
num_orbits = phi_original[-1] / (2 * np.pi)
print('The star does', num_orbits, 'orbits.')
for window_count in range(1, int(num_orbits)):
    try:
        result = minimize(opt, 0.005, method='Powell', args=(window_count,))
        if result.success:
            # print(result.message)
            print('For orbit:', window_count, 'Angle of precession:',
                  result.x)
            # print()
        else:
            print('An error occured:')
            print(result.message)
            print()
            exit(1)
    except ValueError as e:
        print('An exception was caught:')
        print(e)
        exit(1)
print() # to compress the output

