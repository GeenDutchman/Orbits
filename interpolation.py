#!/usr/bin/python3
import numpy as np
from scipy.interpolate import interp1d

#for i in range(len(data)):
    #for j in range(len(data[i])):
        #print(data[i][j], end=' ')
    #print()

for i in range(len(data)):
    if data[i][-1] >= np.pi and data[i][-1] <= 3*np.pi:
        new_row = [data[i][-1], data[i][-2]]
        extractor.append(new_row)
        print(new_row)

extracted = np.array(extractor)

np.savetxt("binary2.dat", extracted)

def copy_piece(theta_column, r_column, start_cut, end_cut):
    theta = []
    r = []
    for row in range(len(theta_column)):
        if theta_column[row] >= start_cut and theta_column[row] <= end_cut:
            theta.append(theta_column[row])
            r.append(r_column[row])
    return theta, r      

def extractor(data):
    theta_original = []
    r_original = []
    for row in range(len(data)):
        theta_original.append(data[row][-1])
        r_original.append(data[row][-2])
    return theta_original, r_original
        

def main(epsilon=0.0):
    data = np.loadtxt('binary1.dat', dtype=np.float64)
    theta_original, r_original = extractor(data)

if __name__ == "__main__":
    # main(sys.argv[1:])
    if len(sys.argv) > 1:
        main(sys.argv[1])
    else:
        main()
