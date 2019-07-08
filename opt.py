#!/usr/bin/python3
import sys
import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import minimize

def find_residual(phi_old, r_old, phi_new, r_new):
    residuals = []
    func = interp1d(phi_new, r_new, kind='cubic')
    for i in range(len(phi_old)):
        r_interp = func(phi_old[i])    
        residuals.append(r_old[i] - r_interp)
    return np.linalg.norm(residuals)
        
def shifter(phi_column, epsilon):
    phi_shifted = [x - (2*np.pi + epsilon) for x in phi_column]
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

def extractor(data):
    phi_original = []
    r_original = []
    for row in range(len(data)):
        phi_original.append(data[row][-2]) # this will now only work for GR
        r_original.append(data[row][-3])
    return phi_original, r_original
        
data = np.loadtxt('binary1.dat', dtype=np.float64)
phi_original, r_original = extractor(data)
phi_old, r_old = copy_piece(phi_original, r_original, np.pi, 3*np.pi)

def opt(epsilon):
    epsilon = epsilon[0]
    phi_shifted = shifter(phi_original, epsilon)
    phi_new, r_new = copy_piece(phi_shifted, r_original, (np.pi/2), ((7*np.pi)/2))
    if phi_old[-1] >= phi_new[-1]:
        raise ValueError('The data needs to have more orbits!!')
    return find_residual(phi_old, r_old, phi_new, r_new)

try:
    num_orbits = phi_original[-1] / (2 * np.pi)
    print('The star does', num_orbits, 'orbits.')
    result = minimize(opt, 0.005, method='Powell')
    if result.success:
        print(result.message)
        print('Angle of precession:', result.x)
        print()
    else:
        print('An error occured:')
        print(result.message)
        print()
        exit(1)
except ValueError as e:
    print('An exception was caught:')
    print(e)
    exit(1)

