"""
  PN EOM Based on Blanchet, "Gravitational Radiation from Post-Newtonian
  Sources and Inspiralling Compact Binaries", Living Rev. Relativity,
  17, (2014), https://link.springer.com/article/10.12942/lrr-2014-2

  N-body PN acceleration from
  Gravity: Newtonian, Post-Newtonian, Relativistic by Poisson and Will
"""

import numpy as np

from RK import RK45_Step

PACKING_ORDER = {
    0: 'Star_pos_x',
    1: 'Star_pos_y',
    2: 'Star_pos_z',
    3: 'Star_vel_x',
    4: 'Star_vel_y',
    5: 'Star_vel_z',
    6: 'Orbit_phase',
    7: 'Binary_separation',
    8: 'Binary_phase',
    9: 'Binary_frequency',
}

PACKING_ORDER_LOOKUP = {
    'Star_pos_x': 0,
    'Star_pos_y': 1,
    'Star_pos_z': 2,
    'Star_vel_x': 3,
    'Star_vel_y': 4,
    'Star_vel_z': 5,
    'Orbit_phase': 6,
    'Binary_separation': 7,
    'Binary_phase': 8,
    'Binary_frequency': 9,
}

N_EQNS = 10

assert(len(PACKING_ORDER) == N_EQNS)
for i in range(N_EQNS):
    assert(PACKING_ORDER_LOOKUP[PACKING_ORDER[i]] == i)


def Omega_of_r(r, **kwargs):
    if 'G' in kwargs:
        G = kwargs['G']
    else:
        G = 1.0

    if 'c_light' in kwargs:
        c_light = kwargs['c_light']
    else:
        c_light = 1.0

    if 'mass' in kwargs:
        mass = kwargs['mass']
    else:
        mass = 1.0

    if 'massratio' in kwargs:
        q = kwargs['massratio']
    else:
        q = 1.0

    nu = q / (1.0 + q)**2  # See Eq. (215). Note, q = m1/m2

    gamma = G * mass / r / c_light**2  # Eq. (225)

    # Omega**2 given by Eq. (228)
    Omega2 = G * mass / r**3 * (
        1.0 + (nu - 3.0) * gamma
            + (6.0 + 41.0/4.0 * nu + nu**2) * gamma**2
            + (-10 + (-75707./840. + 41.0/64. * np.pi**2) * nu
               + 19.0/2.0 * nu**2 + nu**3) * gamma**3
    )

    return np.sqrt(Omega2)


def PN_Acceleration(Xstar, Vstar, rBH, psiBH, OmegaBH, **kwargs):
    if 'G' in kwargs:
        G = kwargs['G']
    else:
        G = 1.0

    if 'c_light' in kwargs:
        c_light = kwargs['c_light']
    else:
        c_light = 1.0

    if 'mass' in kwargs:
        mass = kwargs['mass']
    else:
        mass = 1.0

    if 'massratio' in kwargs:
        q = kwargs['massratio']
    else:
        q = 1.0

    if 'stellar mass' in kwargs:
        Mstar = kwargs['stellar mass']
    else:
        Mstar = 0.0

    M1 = mass * q / (1.0 + q)
    M2 = mass * 1 / (1.0 + q)

    Q1 = M1 / mass
    Q2 = M2 / mass

    Delta = Q1 - Q2

    cospsi = np.cos(psiBH)
    sinpsi = np.sin(psiBH)

    nu = q / (1.0 + q)**2  # See Eq. (215). Note, q = m1/m2

    gamma = G * mass / rBH / c_light**2  # Eq. (225)

    # rdot is given by Eq. (227a)
    rdot = -64.0 / 5.0 * G**3 * mass**3 * nu / rBH**3 / c_light**5 * (
        1.0 + gamma * (-1751.0/336. - 7.0/4.0 * nu)
    )

    # Omega dot is given by Eq. (227b)
    Omegadot = 96.0/5.0 * G * mass * nu / rBH**3 * gamma**(5.0/2.) * (
        1.0 + gamma * (-2591.0/336.0 - 11.0/12.0 * nu)
    )

    # gammadot = d gamma / dr * rdot
    gammadot = -G * mass * rdot / (c_light**2 * rBH**2)

    # rddot = d rdot / dr * rdot
    rddot = (
        7004 * G**3 * gammadot * mass**3 * nu / (105. * c_light**5 * rBH**3)
        + 112 * G**3 * gammadot * mass**3*nu**2 / (5. * c_light**5 * rBH**3)
        + 192 * G**3 * mass**3 * nu * rdot / (5. * c_light**5 * rBH**4)
        - 7004 * G**3 * gamma * mass**3 * nu *
        rdot / (35. * c_light**5 * rBH**4)
        - 336 * G**3 * gamma * mass**3 * nu**2 *
        rdot / (5. * c_light**5 * rBH**4)
    )

    #  Eq.(224a)
    X1 = np.array(
        [
            cospsi*(3*Delta*gamma**2*nu + Q2)*rBH,
            (3*Delta*gamma**2*nu + Q2)*rBH*sinpsi,
            0.0
        ]
    )

    #  Eq.(224b)
    X2 = np.array(
        [
            cospsi*(3*Delta*gamma**2*nu - Q1)*rBH,
            (3*Delta*gamma**2*nu - Q1)*rBH*sinpsi,
            0.0
        ]
    )

    if kwargs['turn_off_gr']:
        X1 = np.array([cospsi*Q2*rBH, Q2*rBH*sinpsi, 0.0])
        X2 = -np.array([cospsi*Q1*rBH, Q1*rBH*sinpsi, 0.0])

    # V1 = d X1(r(t), psi(t)) / dt
    V1 = np.array(
        [
            6*cospsi*Delta*gamma*gammadot*nu*rBH
            + cospsi*(3*Delta*gamma**2*nu + Q2)*rdot
            - OmegaBH*(3*Delta*gamma**2*nu + Q2)*rBH*sinpsi,
            cospsi*OmegaBH*(3*Delta*gamma**2*nu + Q2)*rBH
            + 6*Delta*gamma*gammadot*nu*rBH*sinpsi
            + (3*Delta*gamma**2*nu + Q2)*rdot*sinpsi,
            0.0
        ]
    )

    # V2 = d X2(r(t), psi(t)) / dt
    V2 = np.array(
        [
            6*cospsi*Delta*gamma*gammadot*nu*rBH
            + cospsi*(3*Delta*gamma**2*nu - Q1)*rdot
            - OmegaBH*(3*Delta*gamma**2*nu - Q1)*rBH*sinpsi,
            cospsi*OmegaBH*(3*Delta*gamma**2*nu - Q1)*rBH
            + 6*Delta*gamma*gammadot*nu*rBH*sinpsi
            + (3*Delta*gamma**2*nu - Q1)*rdot*sinpsi,
            0.0
        ]
    )

    # A1 = d V1 / dt
    A1 = np.array(
        [
            6*cospsi*Delta*gammadot**2*nu*rBH
            - cospsi*OmegaBH**2*(3*Delta*gamma**2*nu + Q2)*rBH
            + cospsi*(3*Delta*gamma**2*nu + Q2)*rddot
            + 12*cospsi*Delta*gamma*gammadot*nu*rdot
            - 12*Delta*gamma*gammadot*nu*OmegaBH*rBH*sinpsi
            - 2*OmegaBH*(3*Delta*gamma**2*nu + Q2)*rdot*sinpsi
            - (3*Delta*gamma**2*nu + Q2)*rBH*sinpsi*Omegadot,
            12*cospsi*Delta*gamma*gammadot*nu*OmegaBH*rBH
            + 2*cospsi*OmegaBH*(3*Delta*gamma**2*nu + Q2)*rdot
            + 6*Delta*gammadot**2*nu*rBH*sinpsi
            - OmegaBH**2*(3*Delta*gamma**2*nu + Q2)*rBH*sinpsi
            + (3*Delta*gamma**2*nu + Q2)*rddot*sinpsi
            + 12*Delta*gamma*gammadot*nu*rdot*sinpsi
            + cospsi*(3*Delta*gamma**2*nu + Q2)*rBH*Omegadot,
            0.0
        ]
    )

    # A2 = d V2 / dt
    A2 = np.array(
        [
            6*cospsi*Delta*gammadot**2*nu*rBH
            - cospsi*OmegaBH**2*(3*Delta*gamma**2*nu - Q1)*rBH
            + cospsi*(3*Delta*gamma**2*nu - Q1)*rddot
            + 12*cospsi*Delta*gamma*gammadot*nu*rdot
            - 12*Delta*gamma*gammadot*nu*OmegaBH*rBH*sinpsi
            - 2*OmegaBH*(3*Delta*gamma**2*nu - Q1)*rdot*sinpsi
            - (3*Delta*gamma**2*nu - Q1)*rBH*sinpsi*Omegadot,
            12*cospsi*Delta*gamma*gammadot*nu*OmegaBH*rBH
            + 2*cospsi*OmegaBH*(3*Delta*gamma**2*nu - Q1)*rdot
            + 6*Delta*gammadot**2*nu*rBH*sinpsi
            - OmegaBH**2*(3*Delta*gamma**2*nu - Q1)*rBH*sinpsi
            + (3*Delta*gamma**2*nu - Q1)*rddot*sinpsi
            + 12*Delta*gamma*gammadot*nu*rdot*sinpsi
            + cospsi*(3*Delta*gamma**2*nu - Q1)*rBH*Omegadot,
            0.0
        ]
    )

    X1s = X1 - Xstar
    X2s = X2 - Xstar
    X12 = X1 - X2

    n1s = X1s / np.linalg.norm(X1s)
    n2s = X2s / np.linalg.norm(X2s)
    n12 = X12 / np.linalg.norm(X12)

    Xs1 = -X1s
    Xs2 = -X2s
    X21 = -X12

    ns1 = -n1s
    ns2 = -n2s
    n21 = -n12

    rs1 = np.linalg.norm(Xs1)
    rs2 = np.linalg.norm(Xs2)

    r12 = np.linalg.norm(X12)
    r21 = r12

    r1s = rs1
    r2s = rs2

    # For the Newtonian velocities, assume a circular orbit at sep = rBH
    V1Newton = Q2 * rBH * OmegaBH * np.array((-sinpsi, cospsi, 0.0))
    V2Newton = -Q1 * rBH * OmegaBH * np.array((-sinpsi, cospsi, 0.0))

    # Newtonian accelerations
    A1Newton = G * M2 / rBH**3 * X21
    A2Newton = G * M1 / rBH**3 * X12

    if 'use_newtonian_acceleration_in_PN' in kwargs:
        if kwargs['use_newtonian_acceleration_in_PN']:
            A1 = A1Newton
            A2 = A2Newton

    if 'use_newtonian_velocities_in_PN' in kwargs:
        if kwargs['use_newtonian_velocities_in_PN']:
            V1 = V1Newton
            V2 = V2Newton

    # see https://en.wikipedia.org/wiki/Einstein%E2%80%93Infeld%E2%80%93Hoffmann_equations
    # see Eq. (9.127) in Poisson and Will
    aNewton = (G * M1 * n1s / rs1**2
               + G * M2 * n2s / rs2**2)
    aPN1 = (
        G * M1 * n1s / rs1**2 * (
            np.dot(Vstar, Vstar)
            + 2.0 * np.dot(V1, V1)
            - 4.0 * np.dot(Vstar, V1)
            - 1.5 * np.dot(ns1, V1)**2
            - 4 * G * M1 / rs1
            - 4 * G * M2 / rs2
            - G * M2 / r12
            - 5 * G * Mstar / r1s
            + 0.5 * np.dot(X1s, A1)
        )
        + G * M2 * n2s / rs2**2 * (
            np.dot(Vstar, Vstar)
            + 2.0 * np.dot(V2, V2)
            - 4.0 * np.dot(Vstar, V2)
            - 1.5 * np.dot(ns2, V2)**2
            - 4 * G * M1 / rs1
            - 4 * G * M2 / rs2
            - G * M1 / r21
            - 5 * G * Mstar / r2s
            + 0.5 * np.dot(X2s, A2)
        )
    )

    aPN2 = (
        G * M1 / rs1**2 * (
            np.dot(ns1, 4 * Vstar - 3.0 * V1)*(Vstar - V1)
        )
        +
        G * M2 / rs2**2 * (
            np.dot(ns2, 4 * Vstar - 3.0 * V2)*(Vstar - V2)
        )
    )

    aPN3 = 7.0 / 2.0 * (
        G * M1 * A1 / rs1
        + G * M2 * A2 / rs2
    )

    if kwargs['turn_off_gr']:
        rdot = 0.0
        Omegadot = 0.0
        aPN1 = aPN2 = aPN3 = 0.0
    #  aNewton = -G * mass * Xstar / np.linalg.norm(Xstar)**3

    return aNewton + (aPN1 + aPN2 + aPN3)/c_light**2, rdot, Omegadot


def system_rhs(t, y, **kwargs):

    # unpack
    X_star = y[PACKING_ORDER_LOOKUP['Star_pos_x']               : PACKING_ORDER_LOOKUP['Star_pos_z']+1]
    V_star = y[PACKING_ORDER_LOOKUP['Star_vel_x']               : PACKING_ORDER_LOOKUP['Star_vel_z']+1]
    phi_star = y[PACKING_ORDER_LOOKUP['Orbit_phase']]
    r_BBH = y[PACKING_ORDER_LOOKUP['Binary_separation']]
    phi_BBH = y[PACKING_ORDER_LOOKUP['Binary_phase']]
    Omega_BBH = y[PACKING_ORDER_LOOKUP['Binary_frequency']]

    V_star_dot, r_BBH_dot, Omega_BBH_dot = PN_Acceleration(
        X_star, V_star, r_BBH, phi_BBH, Omega_BBH, **kwargs
    )

    X_star_dot = V_star

    #  omega_star = |L|/r^2 / m_star

    phi_star_dot = np.linalg.norm(
        np.cross(X_star, V_star)) / np.dot(X_star, X_star)

    phi_BBH_dot = Omega_BBH

    ydot = np.zeros(N_EQNS)
    ydot[PACKING_ORDER_LOOKUP['Star_pos_x']         : PACKING_ORDER_LOOKUP['Star_pos_z']+1] = X_star_dot
    ydot[PACKING_ORDER_LOOKUP['Star_vel_x']         : PACKING_ORDER_LOOKUP['Star_vel_z']+1] = V_star_dot
    ydot[PACKING_ORDER_LOOKUP['Orbit_phase']] = phi_star_dot
    ydot[PACKING_ORDER_LOOKUP['Binary_separation']] = r_BBH_dot
    ydot[PACKING_ORDER_LOOKUP['Binary_phase']] = Omega_BBH
    ydot[PACKING_ORDER_LOOKUP['Binary_frequency']] = Omega_BBH_dot

    return ydot


def initial_data(Rstar, rBBH, **kwargs):
    if 'G' in kwargs:
        G = kwargs['G']
    else:
        G = 1.0

    if 'mass' in kwargs:
        mass = kwargs['mass']
    else:
        mass = 1.0

    Sstar = np.sqrt(G * mass / Rstar) * 1.01

    y = np.zeros(N_EQNS)

    Xstar = np.array((Rstar, 0.0, 0.0))
    Vstar = np.array((0.0, Sstar, 0.0))

    y[PACKING_ORDER_LOOKUP['Star_pos_x']        : PACKING_ORDER_LOOKUP['Star_pos_z']+1] = Xstar
    y[PACKING_ORDER_LOOKUP['Star_vel_x']        : PACKING_ORDER_LOOKUP['Star_vel_z']+1] = Vstar
    y[PACKING_ORDER_LOOKUP['Orbit_phase']] = 0
    y[PACKING_ORDER_LOOKUP['Binary_separation']] = rBBH
    y[PACKING_ORDER_LOOKUP['Binary_phase']] = 0.0
    y[PACKING_ORDER_LOOKUP['Binary_frequency']] = Omega_of_r(rBBH)

    return y


def center_of_mass_coordinates_to_BH_positions(r, psi, **kwargs):
    if 'G' in kwargs:
        G = kwargs['G']
    else:
        G = 1.0

    if 'c_light' in kwargs:
        c_light = kwargs['c_light']
    else:
        c_light = 1.0

    if 'mass' in kwargs:
        mass = kwargs['mass']
    else:
        mass = 1.0

    if 'massratio' in kwargs:
        q = kwargs['massratio']
    else:
        q = 1.0

    nu = q / (1.0 + q)**2  # See Eq. (215). Note, q = m1/m2

    gamma = G * mass / r / c_light**2  # Eq. (225)

    Q1 = q / (1.0 + q)  # Q1 = m1/ m
    Q2 = 1.0 / (1.0 + q)  # Q2 = m2/m
    Delta = Q1 - Q2

    rvec = r * np.array((np.cos(psi), np.sin(psi), 0.0))
    y1 = rvec * (Q2 + 3 * gamma**2 * nu * Delta)  # Eq.(224a)
    y2 = rvec * (-Q1 + 3 * gamma**2 * nu * Delta)  # Eq.(224b)

    return y1, y2


def main():
    M = 1.0
    G = 1.0
    C = 1.0
    Q = 1.0

    kwargs = {'G': G, 'mass': M, 'c_light': C, 'massratio': Q,
              'use_newtonian_acceleration_in_PN': True,
              'use_newtonian_velocities_in_PN': False, 'turn_off_gr': False}

    rBBH = 1.0e8
    Rstar = 10.0 * rBBH
    y = initial_data(Rstar, rBBH, **kwargs)

    t = 0.0
    dt = 1.0

    print('# time star_x star_y star_angle star_r')

    tmax = 1.0e15
    while t < tmax:
        t, y, dt = RK45_Step(t, y, dt, system_rhs, **kwargs)
        star_r = np.linalg.norm(y[PACKING_ORDER_LOOKUP['Star_pos_x']:PACKING_ORDER_LOOKUP['Star_pos_z'] + 1])
        print(t, y[0], y[1], y[PACKING_ORDER_LOOKUP['Orbit_phase']], star_r)


if __name__ == "__main__":
    main()
