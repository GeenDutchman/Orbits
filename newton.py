from __future__ import print_function
import numpy as np
from RK import RK4_Step


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

    vdot = -mass * G / r3 * x_vec  # Newton's law of gravitation

    xdot = v_vec

    return np.concatenate((xdot, vdot))

# The Initial condition for the star orbiting a black hole
# The black hole's mass is given by the 'mass',
# Newton's constant by 'G', the initial position by
# (x0, y0, z0), the initial velocity by (vx, vy, vz)


kwargs = {'mass': 1.0, 'G': 1.0}
x0 = 1.0
y0 = 0.0
z0 = 0.0

vx0 = 0.0
vy0 = 1.0
vz0 = 0.0

t = 0.0

# dt is the timestep. The error will be proportional to dt**4
dt = 1.0e-3

initial_position = np.array((x0, y0, z0), dtype=np.float64)
initial_velocity = np.array((vx0, vy0, vz0), dtype=np.float64)

Y = np.concatenate((initial_position, initial_velocity))


while t < 10:
    print(Y[0], Y[1], Y[2])

    # The Runge-Kutta routine returns the new value of Y, t, and a
    # possibly updated value of dt
    t, Y, dt = RK4_Step(t, Y, dt, Keppler_RHS, **kwargs)
