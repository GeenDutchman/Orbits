"""
  PN EOM Based on Blanchet, "Gravitational Radiation from Post-Newtonian
  Sources and Inspiralling Compact Binaries", Living Rev. Relativity,
  17, (2014),

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

  nu = q/(1.0+q)**2 # See Eq. (215). Note, q = m1/m2

  gamma = G * mass / r / c_light**2 # Eq. (225)

  # Omega**2 given by Eq. (228)
  Omega2 = (G * mass / r**3) * ( 1.0 + (nu - 3.0) * gamma + (6.0
    +41.0/4.0 * nu + nu**2) * gamma**2 + (-10 + (-75707./840. +
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

  nu = q/(1.0+q)**2 # See Eq. (215). Note, q = m1/m2

  gamma = G * mass / r / c_light**2 # Eq. (225)

  # rdot is given by Eq. (227a)
  rdot = -64.0 / 5.0 * G**3 * mass**3 * nu / r**3 / c_light**5 * (1.0
      + gamma * (-1751.0/336. - 7.0/4.0 * nu))

  # Omega dot is given by Eq. (227b)
  Omegadot = 96.0/5.0 * G * mass * nu / r**3 * gamma**(5.0/2.) * (1.0
      + gamma * (-2591.0/336.0 - 11.0/12.0 * nu))

  # take psi_dot = Omega

  psidot = Omega

  return np.array((rdot, psidot, Omegadot))


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

  nu = q/(1.0+q)**2 # See Eq. (215). Note, q = m1/m2

  gamma = G * mass / r / c_light**2 # Eq. (225)


  X1 = q / (1.0 + q) #  X1 = m1/ m
  X2 = 1.0 / (1.0 + q) # X2 = m2/m
  Delta = X1 - X2

  rvec = r * np.array((np.cos(psi), np.sin(psi), 0.0))
  y1 = rvec * (X2 + 3 * gamma**2 * nu * Delta) # Eq.(224a)
  y2 = rvec * (-X1 + 3 * gamma**2 * nu * Delta) # Eq.(224b)

  return y1, y2


if __name__ == "__main__":
  M = 1.0
  G = 1.0
  C = 1.0

  R=100.0

  Q = 1.0

  kwargs={'G' : G, 'mass' : M, 'c_light' : C, 'massratio' : Q,
      'max_step_size' : 100}

  Omega = Omega_of_r(R, **kwargs)

  psi = 0.0

  y = np.array((R, psi, Omega))

  t = 0.0

  dt = 1.0e-1
  while y[0] > 10:
    x1, x2 = center_of_mass_coordinates_to_BH_positions(y[0], y[1], **kwargs)
    print (t, x1[0] )

    t, y, dt = RK.RK45_Step(t, y, dt, PN_Orbit, **kwargs)

