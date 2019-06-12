from  __future__ import print_function
import numpy as np
from RK import RK4_Step, RK45_Step


def Keppler_RHS(t, y0, **kwargs):
  """
     Outputs the time derivatives of the position vector x and
     velocity vector v according to Newton's law of gravitation.
     The inputs are a time t, which is ignored, and a numpy vector y0.
     y0 should contain 6 elements. The first three are the components
     of x, the last three are the components of v. The output is a
     numpy array with six elements corresponding to the time
     derivatives of the components of x and v. The optional kwargs can
     contain mass and Newton's constant parameters.
   """

  x_vec = y0[0:3]  # position stored in first 3 elements
  v_vec = y0[3:]   # velocity stored in last 3 elements

  r = np.linalg.norm(x_vec)
  r3 = r**3

  if 'mass' in kwargs:
    mass = kwargs['mass']
  else:
    mass = 1.0

  if 'G' in kwargs:
    G = kwargs['G']
  else:
    G = 1.0
 
  vdot = -mass * G / r3 * x_vec # Newton's law of gravitation

  xdot = v_vec

  return np.concatenate((xdot, vdot))

### The Initial condition for the star orbiting a black hole
### The black hole's mass is given by the 'mass',
### Newton's constant by 'G', the initial position by
### (x0, y0, z0), the initial velocity by (vx, vy, vz)

kwargs = { 'mass' : 1.0, 'G' : 1.0 }
x0 = 1.0
y0 = 0.0
z0 = 0.0

vx0 = 0.0
vy0 = 0.9
vz0 = 0.0

t = 0.0

# dt is the timestep. The error will be proportional to dt**4
dt = 1.0e-3

tmax = 100

initial_position = np.array((x0, y0, z0), dtype=np.float64)
initial_velocity = np.array((vx0, vy0, vz0), dtype=np.float64)

Y = np.concatenate((initial_position, initial_velocity))

list_size = int(tmax / dt) + 1

r_list = np.zeros(shape=( list_size, ))
phi_list = np.zeros(shape=( list_size, ))
x_list = np.zeros(shape=( list_size, ))
y_list = np.zeros(shape=( list_size, ))
t_list = np.zeros(shape=( list_size, ))

lst_indx = 0
while t < tmax:
  x_pos = Y[0]
  y_pos = Y[1]
  z_pos = Y[2]
  r = np.sqrt(x_pos**2 + y_pos**2)
  phi = np.arctan2(y_pos, x_pos)
  phi_list[lst_indx] = phi
  r_list[lst_indx] = r
  x_list[lst_indx] = Y[0]
  y_list[lst_indx] = Y[1]
  t_list[lst_indx] = t
  lst_indx += 1
  #print (t, Y[0], Y[1], r, phi)

  ## The Runge-Kutta routine returns the new value of Y, t, and a
  ## possibly updated value of dt
  t, Y, dt = RK4_Step(t, Y, dt, Keppler_RHS, **kwargs)

phi =  np.unwrap(phi_list[0:lst_indx])
for i in range(lst_indx):
  print (t_list[i], x_list[i], y_list[i], r_list[i], phi[i])

