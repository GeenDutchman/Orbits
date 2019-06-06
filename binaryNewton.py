from __future__ import print_function
import numpy as np
from RK import RK4_Step, RK45_Step


def Keppler_Binary_RHS(t, y0, BHs0, **kwargs):
    star_x_vec = y0[0:3] # star position stored in first 3 elements
    star_v_vec = y0[3:] # star velocity stored in last 3 elements

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

    BHs_x_init_vec = BHs0[0:3]
    BHs_v_init_vec = BHs0[3:]

    # calculate the current position, but does not do the z coord??
    BHs_x_curr_vec = BHs_x_init_vec * np.array((np.sin(BHs_v_init_vec[0] * t), np.cos(BHs_v_init_vec[1] * t), 0),
                                               dtype=np.float64)  # , np.cos(BHs_v_init_vec[2] * t)))

    BH1_mass = combined_BH_mass / (BH_ratio + 1)
    BH2_mass = combined_BH_mass - BH1_mass  # (-1 * combined_BH_mass * BH_ratio) / (BH_ratio + 1)

    BH1_pos =
