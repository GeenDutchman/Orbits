#!/usr/bin/python3
import numpy as np
from scipy.interpolate import interp1d

data = np.loadtxt('binary1.dat', dtype=np.float64)
#for i in range(len(data)):
    #for j in range(len(data[i])):
        #print(data[i][j], end=' ')
    #print()

extractor = []

for i in range(len(data)):
    if data[i][-1] >= np.pi and data[i][-1] <= 3*np.pi:
        new_row = [data[i][-1], data[i][-2]]
        extractor.append(new_row)
        print(new_row)

extracted = np.array(extractor)

np.savetxt("binary2.dat", extracted)




