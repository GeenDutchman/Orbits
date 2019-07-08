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

def extractor(data):
    phi_original = []
    r_original = []
    for row in range(len(data)):
        phi_original.append(data[row][-2]) # this will now only work for GR
        r_original.append(data[row][-3])
    return phi_original, r_original
        
data = np.loadtxt('binary1.dat', dtype=np.float64)
phi_original, r_original = extractor(data)
phi_old, r_old = copy_piece(phi_original, r_original, 2*np.pi - np.pi, 2*np.pi + np.pi)

def opt(epsilon, *args):
    nOrbits = args[0]
    epsilon = epsilon[0]
    phi_shifted = shifter(phi_original, epsilon, nOrbits)
    phi_new, r_new = copy_piece(phi_shifted, r_original,2*np.pi - 1.5*np.pi,2*np.pi + 1.5*np.pi)
    if phi_old[-1] >= phi_new[-1]:
        raise ValueError('The data needs to have more orbits!!')
    return find_residual(phi_old, r_old, phi_new, r_new)

num_orbits = phi_original[-1] / (2 * np.pi)
print('The star does', num_orbits, 'orbits.')
for nOrbits in range(1, int(num_orbits)):
    try:
        result = minimize(opt, 0.005 , method='Powell', args=(nOrbits,))
        if result.success:
            print(result.message)
            print('For orbit:',nOrbits + 1,'Angle of precession:', result.x)
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

