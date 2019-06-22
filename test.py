import numpy as np
from scipy.interpolate import interp1d

LEN=10
x=np.zeros(LEN)

for i in range(LEN):
  x[i] = 0 + i * 2*np.pi/LEN


y = np.sin(x)

z = np.ndarray(shape=(LEN,2), dtype=np.float64)

for i in range(LEN):
  z[i] = (x[i], y[i])

np.savetxt("outout", z)


LEN2=100
x2 = np.zeros(LEN2)
z2 = np.zeros(shape=(LEN2,2), dtype=np.float64)

f2 = interp1d(x, y, kind='cubic')
for i in range(LEN2-20):
  x2[i] =  (i+10) * 2*np.pi/LEN2
  


for i in range(LEN2):
  z2[i] = (x2[i], f2(x2[i]))

np.savetxt("outout2", z2)

