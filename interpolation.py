#!/usr/bin/python3
import sys
import numpy as np
from scipy.interpolate import interp1d

def interpolater(phi_old, r_old, phi_new, r_new):
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
        if phi_column[row] >= start_cut and phi_column[row] <= end_cut:
            phi.append(phi_column[row])
            r.append(r_column[row])
    return phi, r      

def extractor(data):
    phi_original = []
    r_original = []
    for row in range(len(data)):
        phi_original.append(data[row][-1])
        r_original.append(data[row][-2])
    return phi_original, r_original
        
def main(epsilon=0.0):
    data = np.loadtxt('binary1.dat', dtype=np.float64)
    phi_original, r_original = extractor(data)
    phi_old, r_old = copy_piece(phi_original, r_original, np.pi, 3*np.pi)
    phi_shifted = shifter(phi_original, epsilon)
    phi_new, r_new = copy_piece(phi_shifted, r_original, (np.pi/2), ((7*np.pi)/2))
    residual = interpolater(phi_old, r_old, phi_new, r_new)
    print('Residual:')
    print(residual)

if __name__ == "__main__":
    # main(sys.argv[1:])
    if len(sys.argv) > 1:
        main(float(sys.argv[1]))
    else:
        main()
