"""
  PN EOM Based on Blanchet, "Gravitational Radiation from Post-Newtonian
  Sources and Inspiralling Compact Binaries", Living Rev. Relativity,
  17, (2014), https://link.springer.com/article/10.12942/lrr-2014-2#PartA

"""

import numpy as np
import RK


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

    nu = q/(1.0+q)**2  # See Eq. (215). Note, q = m1/m2

    gamma = G * mass / r / c_light**2  # Eq. (225)

    # Omega**2 given by Eq. (228)
    Omega2 = (G * mass / r**3) * (1.0 + (nu - 3.0) * gamma + (6.0
                                                              + 41.0/4.0 * nu + nu**2) * gamma**2 + (-10 + (-75707./840. +
                                                                                                            41.0/64. * np.pi**2) * nu + 19.0/2.0 * nu**2 + nu**3) * gamma**3)

    return np.sqrt(Omega2)


def PN_Orbit(t, r, psi, Omega, **kwargs):
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

    nu = q/(1.0+q)**2  # See Eq. (215). Note, q = m1/m2

    gamma = G * mass / r / c_light**2  # Eq. (225)

    # rdot is given by Eq. (227a)
    rdot = -64.0 / 5.0 * G**3 * mass**3 * nu / r**3 / c_light**5 * (1.0
                                                                    + gamma * (-1751.0/336. - 7.0/4.0 * nu))

    # Omega dot is given by Eq. (227b)
    Omegadot = 96.0/5.0 * G * mass * nu / r**3 * gamma**(5.0/2.) * (1.0
                                                                    + gamma * (-2591.0/336.0 - 11.0/12.0 * nu))

    # take psi_dot = Omega

    psidot = Omega

    return np.array((rdot, psidot, Omegadot))


def PN_Acceleration(Xstar, Vstar, rBH, psiBH, Omega, **kwargs):
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

    M1 = mass * q / (1.0 + q)
    M2 = mass / (1.0 + q)

    Q1 = M1 / mass
    Q2 = M2 / mass

    Delta = Q1 - Q2

    cospsi = np.cos(psiBH)
    sinpsi = np.sin(psiBH)

    nu = q/(1.0+q)**2  # See Eq. (215). Note, q = m1/m2

    gamma = G * mass / rBH / c_light**2  # Eq. (225)

    # rdot is given by Eq. (227a)
    rdot = -64.0 / 5.0 * G**3 * mass**3 * nu / rBH**3 / c_light**5 * (1.0
                                                                      + gamma * (-1751.0/336. - 7.0/4.0 * nu))

    # Omega dot is given by Eq. (227b)
    Omegadot = 96.0/5.0 * G * mass * nu / rBH**3 * gamma**(5.0/2.) * (1.0
                                                                      + gamma * (-2591.0/336.0 - 11.0/12.0 * nu))

    gammadot = -((G*mass*rdot)/(c_light**2*rBH**2))

    rddot = (7004*G**3*gammadot*mass**3*nu)/(105.*c_light**5*rBH**3) + \
            (112*G**3*gammadot*mass**3*nu**2)/(5.*c_light**5*rBH**3) + \
            (192*G**3*mass**3*nu*rdot)/(5.*c_light**5*rBH**4) - \
            (7004*G**3*gamma*mass**3*nu*rdot)/(35.*c_light**5*rBH**4) - \
            (336*G**3*gamma*mass**3*nu**2*rdot)/(5.*c_light**5*rBH**4)

    X1 = np.array(
        [cospsi*(3*Delta*gamma**2*nu + Q2)*rBH, (3*Delta*gamma**2*nu +
                                                 Q2)*rBH*sinpsi, 0]
    )

    X2 = np.array(
        [cospsi*(3*Delta*gamma**2*nu - Q1)*rBH, (3*Delta*gamma**2*nu -
                                                 Q1)*rBH*sinpsi, 0]
    )

    V1 = np.array(
        [6*cospsi*Delta*gamma*gammadot*nu*rBH + cospsi*(3*Delta*gamma**2*nu + Q2)*rdot -
         Omega*(3*Delta*gamma**2*nu + Q2)*rBH*sinpsi,
         cospsi*Omega*(3*Delta*gamma**2*nu + Q2)*rBH + 6*Delta*gamma*gammadot*nu*rBH*sinpsi +
         (3*Delta*gamma**2*nu + Q2)*rdot*sinpsi, 0]
    )

    V2 = np.array(
        [6*cospsi*Delta*gamma*gammadot*nu*rBH + cospsi*(3*Delta*gamma**2*nu - Q1)*rdot -
         Omega*(3*Delta*gamma**2*nu - Q1)*rBH*sinpsi,
         cospsi*Omega*(3*Delta*gamma**2*nu - Q1)*rBH + 6*Delta*gamma*gammadot*nu*rBH*sinpsi +
         (3*Delta*gamma**2*nu - Q1)*rdot*sinpsi, 0]
    )

    A1 = np.array(
        [6*cospsi*Delta*gammadot**2*nu*rBH - cospsi*Omega**2*(3*Delta*gamma**2*nu + Q2)*rBH +
         cospsi*(3*Delta*gamma**2*nu + Q2)*rddot + 12*cospsi*Delta*gamma*gammadot*nu*rdot -
         12*Delta*gamma*gammadot*nu*Omega*rBH*sinpsi - 2*Omega*(3*Delta*gamma**2*nu + Q2)*rdot*sinpsi -
         (3*Delta*gamma**2*nu + Q2)*rBH*sinpsi*Omegadot,
         12*cospsi*Delta*gamma*gammadot*nu*Omega*rBH + 2*cospsi*Omega*(3*Delta*gamma**2*nu + Q2)*rdot +
         6*Delta*gammadot**2*nu*rBH*sinpsi - Omega**2*(3*Delta*gamma**2*nu + Q2)*rBH*sinpsi +
         (3*Delta*gamma**2*nu + Q2)*rddot*sinpsi + 12*Delta*gamma*gammadot*nu*rdot*sinpsi +
         cospsi*(3*Delta*gamma**2*nu + Q2)*rBH*Omegadot, 0]
    )

    A2 = np.array(
        [6*cospsi*Delta*gammadot**2*nu*rBH - cospsi*Omega**2*(3*Delta*gamma**2*nu - Q1)*rBH +
         cospsi*(3*Delta*gamma**2*nu - Q1)*rddot + 12*cospsi*Delta*gamma*gammadot*nu*rdot -
         12*Delta*gamma*gammadot*nu*Omega*rBH*sinpsi - 2*Omega*(3*Delta*gamma**2*nu - Q1)*rdot*sinpsi -
         (3*Delta*gamma**2*nu - Q1)*rBH*sinpsi*Omegadot,
         12*cospsi*Delta*gamma*gammadot*nu*Omega*rBH + 2*cospsi*Omega*(3*Delta*gamma**2*nu - Q1)*rdot +
         6*Delta*gammadot**2*nu*rBH*sinpsi - Omega**2*(3*Delta*gamma**2*nu - Q1)*rBH*sinpsi +
         (3*Delta*gamma**2*nu - Q1)*rddot*sinpsi + 12*Delta*gamma*gammadot*nu*rdot*sinpsi +
         cospsi*(3*Delta*gamma**2*nu - Q1)*rBH*Omegadot, 0]
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

    V1Newton = Q2 * rBH * Omega * np.array((-sinpsi, cospsi, 0))
    V2Newton = -Q1 * rBH * Omega * np.array((-sinpsi, cospsi,0))

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
    aNewton = G * M1 * n1s / rs1**2 + G * M2 * n2s / rs2**2
    aPN1 = (
        G * M1 * n1s / rs1**2 * (np.dot(Vstar, Vstar) +
                                 2.0 * np.dot(V1, V1) - 4.0 * np.dot(Vstar, V1) - 1.5 * np.dot(ns1,
                                                                                               V1)**2 - 4 * G * M1 / rs1 - 4 * G * M2 / rs2 - G * M2 / r12 +
                                 0.5 * np.dot(X1s, A1))
        +
        G * M2 * n2s / rs2**2 * (np.dot(Vstar, Vstar) +
                                 2.0 * np.dot(V2, V2) - 4.0 * np.dot(Vstar, V2) - 1.5 * np.dot(ns2,
                                                                                               V2)**2 - 4 * G * M1 / rs1 - 4 * G * M2 / rs2 - G * M1 / r12 +
                                 0.5 * np.dot(X2s, A2)))

    aPN2 = (
        G * M1 / rs1**2 * (np.dot(ns1, 4 * Vstar - 3.0 * V1)*(Vstar -
                                                              V1))
        +
        G * M2 / rs2**2 * (np.dot(ns2, 4 * Vstar - 3.0 * V2)*(Vstar -
                                                              V2)))

    aPN3 = 7.0 / 2.0 * (G * M1 * A1 / rs1 + G * M2 * A2 / rs2)

    # Note that this will override use_newtonian_velocities_in_PN and use_newtonian_acceleration_in_PN
    if 'use_post_newtonian_corrections_in_PN' in kwargs:
        if kwargs['use_post_newtonian_corrections_in_PN'] is False:
            aPN1 *= 0
            aPN2 *= 0
            aPN3 *= 0

    return aNewton + (aPN1 + aPN2 + aPN3)/c_light**2, rdot, Omegadot


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

    nu = q/(1.0+q)**2  # See Eq. (215). Note, q = m1/m2

    gamma = G * mass / r / c_light**2  # Eq. (225)

    X1 = q / (1.0 + q)  # X1 = m1/ m
    X2 = 1.0 / (1.0 + q)  # X2 = m2/m
    Delta = X1 - X2

    rvec = r * np.array((np.cos(psi), np.sin(psi), 0.0))
    y1 = rvec * (X2 + 3 * gamma**2 * nu * Delta)  # Eq.(224a)
    y2 = rvec * (-X1 + 3 * gamma**2 * nu * Delta)  # Eq.(224b)

    return y1, y2


if __name__ == "__main__":
    M = 1.0
    G = 1.0
    C = 1.0

    R = 100.0

    Q = 1.0

    kwargs = {'G': G, 'mass': M, 'c_light': C, 'massratio': Q,
              'max_step_size': 100}

    Omega = Omega_of_r(R, **kwargs)

    psi = 0.0

    y = np.array((R, psi, Omega))

    t = 0.0

    dt = 1.0e-1
    while y[0] > 10:
        x1, x2 = center_of_mass_coordinates_to_BH_positions(
            y[0], y[1], **kwargs)
        print(t, x1[0])

        t, y, dt = RK.RK45_Step(t, y, dt, PN_Orbit, **kwargs)
