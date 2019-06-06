from __future__ import print_function
import numpy as np
from RK import RK4_Step, RK45_Step


def Keppler_Binary_RHS(t, y0, BHs0, **kwargs):
    star_x_vec = y0[0:3]  # star position stored in first 3 elements
    star_v_vec = y0[3:]  # star velocity stored in last 3 elements

    # star position vector norm, such that
    # star_r = (star_x**2 + star_y**2 + star_z**2)**1/2
    # or the distance from the origin
    star_r = np.linalg.norm(star_x_vec)

    if 'mass' in kwargs:
        combined_BH_mass = kwargs['mass']
    else:
        combined_BH_mass = 1.0

    if 'G' in kwargs:
        G = kwargs['G']
    else:
        G = 1.0

    if 'q' in kwargs:
        BH_ratio = kwargs['q']
    else:
        BH_ratio = 1.0  # equal mass default

    if 'bhs' in kwargs:
        BHs0 = kwargs['bhs']
    else:
        print("Must have black hole data!!")
        exit(2)

    BHs_x_init_vec = BHs0[0:3]
    BHs_v_init_vec = BHs0[3:]

    # calculate the current position, but does not do the z coord??
    BHs_x_curr_vec = BHs_x_init_vec * np.array((np.sin(BHs_v_init_vec[0] * t), np.cos(BHs_v_init_vec[1] * t), 0),
                                               dtype=np.float64)  # , np.cos(BHs_v_init_vec[2] * t)))

    BH1_mass = combined_BH_mass / (BH_ratio + 1)
    # (-1 * combined_BH_mass * BH_ratio) / (BH_ratio + 1)
    BH2_mass = combined_BH_mass - BH1_mass

    BH1_pos = BHs_x_curr_vec / (BH_ratio + 1)
    BH2_pos = BH1_pos - BHs_x_curr_vec

    # find the acceleration due to gravity from the respective blach holes
    acc_star_1 = (star_r - BH1_pos) * (-1 * BH1_mass * G) / \
        (np.linalg.norm(star_r - BH1_pos) ** 3)
    acc_star_2 = (star_r - BH2_pos) * (-1 * BH2_mass * G) / \
        (np.linalg.norm(star_r - BH2_pos) ** 3)

    acc_star = acc_star_1 + acc_star_2
    vel_star = star_v_vec

    return np.concatenate((vel_star, acc_star))

# The Initial condition for the star orbiting a black hole
# The black hole's mass is given by the 'mass',
# Newton's constant by 'G', the initial position by
# (x0, y0, z0), the initial velocity by (vx, vy, vz)


kwargs = {'mass': 1.0, 'G': 1.0, 'q': 1.0}
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

BHx = 0.0
BHy = 2.0
BHz = 0.0

BHvx = 1.0
BHvy = 0.0
BHvz = 0.0

initial_bh_pos = np.array((BHx, BHy, BHz), dtype=np.float64)
initial_bh_velocity = np.array((BHvx, BHvy, BHvz), dtype=np.float64)

BHs0 = np.concatenate((initial_bh_pos, initial_bh_velocity))

kwargs['bhs'] = BHs0


while t < 30:
    print(Y[0], Y[1], Y[2])

    # The Runge-Kutta routine returns the new value of Y, t, and a
    # possibly updated value of dt
    t, Y, dt = RK4_Step(t, Y, dt, Keppler_Binary_RHS, **kwargs)
