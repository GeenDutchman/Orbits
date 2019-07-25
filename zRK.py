import numpy as np
""" 
A very simple RK4 and RK45 implementation 
"""


DEFAULT_TOL = 1.0e-6
DEFAULT_MAX_STEP_SIZE = 1.0e100

def RK4_Step(xold, yold, dx, fxy, **kwargs):
  """
  Implement NON-adaptive RK4
  Arguments:
    xold : Old value of independent variable. Should be a scalar value
           of type float, or np.float64.

    yold : Old value (vector) of dependent variables in vector form.
           Should be an np.array of type np.float64 or a scalar of
           type float of np.float64.

    dx : requested stepsize. Should be a scalar of type float or
         np.float64

    fxy: The RHS function. This function should accept arguments of
        type xold, yold, and, optionally, **kwargs

    **kwargs : a dictionary of optional parameters for fxy

  Returns:
    Updated values of x (independent variable), y (vector), dx
    (always the same as the input dx)
  """

  k1 = fxy(xold, yold, **kwargs)
  k2 = fxy(xold+ 0.5 * dx, yold + 0.5 * dx * k1, **kwargs)
  k3 = fxy(xold+ 0.5 * dx, yold + 0.5 * dx * k2, **kwargs)
  k4 = fxy(xold+ 1.0 * dx, yold + 1.0 * dx * k3, **kwargs)

  xnew = xold + dx
  ynew = yold + (dx / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
  return xnew, ynew, dx



# Adaptive RK45
def RK45_Step(xold, yold, dx, fxy, tol=DEFAULT_TOL,
    max_step_size=DEFAULT_MAX_STEP_SIZE, **kwargs):
  """
  Implement Adaptive RK45 (Runge-Kutta-Fehlberg)
  Arguments:
    xold : Old value of independent variable. Should be a scalar value
           of type float, or np.float64.

    yold : Old value (vector) of dependent variables in vector form.
           Should be an np.array of type np.float64 or a scalar of
           type float of np.float64.

    dx : requested stepsize. Should be a scalar of type float or
         np.float64

    fxy: The RHS function. This function should accept arguments of
        type xold, yold, and, optionally, **kwargs

    tol : error tolerance. This is optional. Should be of
          type float of np.float64

    **kwargs : a dictionary of optional parameters for fxy

  Returns:
    Updated values of x (independent variable), y (vector), dx
    (optimal step size)

    If the timestep is rejected because dx is too large, then the
    return values are xold, yold, and dxnew. Where dxnew is an
    approximation to the optimal step size that should be used.
  """


  k1 = dx * fxy(xold, yold, **kwargs)
  k2 = dx * fxy(xold + 1.0/4.0 * dx, yold + 1.0/4.0 * k1, **kwargs)
  k3 = dx * fxy(xold + 3.0/8.0 * dx, yold + 3.0/32.0 * k1
               + 9.0/32.0 * k2, **kwargs)
  k4 = dx * fxy(xold + 12.0/13.0 * dx, yold + 1932.0/2197.0 * k1
               - 7200.0/2197.0 * k2 + 7296.0/2197.0 * k3, **kwargs)
  k5 = dx * fxy(xold + dx, yold + 439.0/216.0 * k1 - 8.0 * k2 
               + 3680.0/513.0 * k3 - 845.0/4104.0 * k4, **kwargs)
  k6 = dx * fxy(xold + 1.0/2.0 * dx, yold -8.0/27.0 * k1
               + 2.0 * k2 - 3544.0/2565.0 * k3 + 1859.0/4104.0 * k4
               - 11.0/40.0 * k5, **kwargs)

  # higher order accuracy ynew
  ynew = yold + 16.0/135.0 * k1 + 6656.0/12825.0 * k3 +\
             28561.0/56430.0 * k4 -9.0/50.0 * k5 + 2.0/55.0 * k6

  # lower order accuracy ynew
  ynew4 = yold + 25.0/216.0 * k1 + 1408.0/2565.0 * k3 +\
             2197.0/4104.0 * k4 -1.0/5.0 * k5

  xnew = xold + dx

  err= abs(np.array(ynew-ynew4)).max()

  if err < tol * 1.0e-3:
    dxnew = 2* dx
  else:
    dxnew = 0.9 * dx * (tol/err)**(0.25)

  if dxnew > max_step_size:
    dxnew = max_step_size 

  # If the error is too big, we need to reject the step.
  # Pass the old values of x and y, but a new smaller value of dx,
  # to the calling function.
  if err > tol:
    return xold, yold, dxnew

  return (xnew, ynew, dxnew)
