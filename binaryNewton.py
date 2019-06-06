from __future__ import print_function
import numpy as np
from RK import RK4_Step, RK45_Step


def calc_omega(mass, G, pos1, pos2):
    magnitude1 = np.linalg.norm(pos1)
    magnitude2 = np.linalg.norm(pos2)
    mag_sum = magnitude1 + magnitude2
    dist = abs(mag_sum)
    dist3 = dist ** 3

    omega = np.sqrt((mass * G) / dist3)
    return omega, magnitude1, magnitude2


def Keppler_Binary_RHS(t, y0, **kwargs):
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

    if 'bh1' in kwargs:
        BH1 = kwargs['bh1']
    else:
        print("Must have black hole data!!")
        exit(2)

    if 'bh2' in kwargs:
        BH2 = kwargs['bh2']
    else:
        print("Must have black hole information!!")
        exit(2)

    BH1_x_vec = BH1[0:3]
    #BH1_v_vec = BH1[3:]

    BH2_x_vec = BH2[0:3]
    #BH2_v_vec = BH2[3:]

    omega, BH1r, BH2r = calc_omega(combined_BH_mass, G, BH1_x_vec, BH2_x_vec)

    # calculate the current position, but does not do the z coord??
    BH1_x_vec[0] = BH1r * np.cos(omega * t)
    BH1_x_vec[1] = BH1r * np.sin(omega * t)
    BH1_x_vec[2] = 0  # don't do z...

    # calculate the current position, but does not do the z coord??
    BH2_x_vec[0] = -1 * BH2r * np.cos(omega * t)
    BH2_x_vec[1] = -1 * BH2r * np.sin(omega * t)
    BH2_x_vec[2] = 0  # don't do z...

    BH1_mass = combined_BH_mass / (BH_ratio + 1)
    # (-1 * combined_BH_mass * BH_ratio) / (BH_ratio + 1)
    BH2_mass = combined_BH_mass - BH1_mass

    # find the acceleration due to gravity from the respective black holes
    acc_star_1 = (star_r - BH1_x_vec) * (-1 * BH1_mass * G) / \
        (np.linalg.norm(star_r - BH1_x_vec) ** 3)
    acc_star_2 = (star_r - BH2_x_vec) * (-1 * BH2_mass * G) / \
        (np.linalg.norm(star_r - BH2_x_vec) ** 3)

    acc_star = acc_star_1 + acc_star_2
    vel_star = star_v_vec

    kwargs['bh1'] = BH1_x_vec
    kwargs['bh2'] = BH2_x_vec

    return np.concatenate((vel_star, acc_star))

# The Initial condition for the star orbiting a black hole
# The black hole's mass is given by the 'mass',
# Newton's constant by 'G', the initial position by
# (x0, y0, z0), the initial velocity by (vx, vy, vz)


kwargs = {'mass': 1.0, 'G': 1.0, 'q': 1.0}
x0 = 2.0
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

BH1x = 1.0
BH1y = 0.0
BH1z = 0.0

'''
BH1vx = 0.0
BH1vy = 1.0
BH1vz = 0.0
'''

initial_bh1_pos = np.array((BH1x, BH1y, BH1z), dtype=np.float64)
#initial_bh1_velocity = np.array((BH1vx, BH1vy, BH1vz), dtype=np.float64)

BH1 = initial_bh1_pos # np.concatenate((initial_bh1_pos, initial_bh1_velocity))

BH2x = -1.0
BH2y = 0.0
BH2z = 0.0

'''
BH2vx = 0.0
BH2vy = -1.0
BH2vz = 0.0
'''

initial_bh2_pos = np.array((BH2x, BH2y, BH2z), dtype=np.float64)
#initial_bh2_velocity = np.array((BH2vx, BH2vy, BH2vz), dtype=np.float64)

BH2 = initial_bh2_pos # np.concatenate((initial_bh2_pos, initial_bh2_velocity))

kwargs['bh1'] = BH1
kwargs['bh2'] = BH2

while t < 30:
    BH1 = kwargs['bh1']
    BH2 = kwargs['bh2']
    print(Y[1], Y[2], BH1[0], BH1[1], BH2[0], BH2[1])

    # The Runge-Kutta routine returns the new value of Y, t, and a
    # possibly updated value of dt
    t, Y, dt = RK4_Step(t, Y, dt, Keppler_Binary_RHS, **kwargs)
