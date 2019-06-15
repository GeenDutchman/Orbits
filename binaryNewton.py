from __future__ import print_function
import numpy as np
from RK import RK4_Step, RK45_Step
import sys
import getopt


def calc_omega(mass, G, pos1, pos2):
    magnitude1 = np.linalg.norm(pos1)
    magnitude2 = np.linalg.norm(pos2)
    mag_sum = magnitude1 + magnitude2
    dist = abs(mag_sum)
    dist3 = dist ** 3

    omega = np.sqrt((mass * G) / dist3)
    return omega, dist


def Keppler_Binary_RHS(t, y0, **kwargs):
    star_x_vec = y0[0:3]  # star position stored in first 3 elements
    star_v_vec = y0[3:]  # star velocity stored in last 3 elements

    # star position vector norm, such that
    # star_r = (star_x**2 + star_y**2 + star_z**2)**1/2
    # or the distance from the origin
    # star_r = np.linalg.norm(star_x_vec)

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

    BH2_x_vec = BH2[0:3]

    if 'omega' in kwargs and 'BH_dist' in kwargs:
        omega = kwargs['omega']
        BH_dist = kwargs['BH_dist']
    else:
        print('# Calculating omega')
        omega, BH_dist = calc_omega(combined_BH_mass, G, BH1_x_vec, BH2_x_vec)
        kwargs['omega'] = omega
        kwargs['BH_dist'] = BH_dist

    half_BH_dist = 0.5 * BH_dist

    # calculate the current position, but does not do the z coord??
    BH1_x_vec[0] = half_BH_dist * np.cos(omega * t)
    BH1_x_vec[1] = half_BH_dist * np.sin(omega * t)
    BH1_x_vec[2] = 0  # don't do z...

    # calculate the current position, but does not do the z coord??
    BH2_x_vec[0] = -1 * half_BH_dist * np.cos(omega * t)
    BH2_x_vec[1] = -1 * half_BH_dist * np.sin(omega * t)
    BH2_x_vec[2] = 0  # don't do z...

    BH1_mass = combined_BH_mass / (BH_ratio + 1)
    # (-1 * combined_BH_mass * BH_ratio) / (BH_ratio + 1)
    BH2_mass = combined_BH_mass - BH1_mass

    # find the acceleration due to gravity from the respective black holes
    acc_star_1 = (star_x_vec - BH1_x_vec) * (-1 * BH1_mass * G) / \
        (np.linalg.norm(star_x_vec - BH1_x_vec) ** 3)
    acc_star_2 = (star_x_vec - BH2_x_vec) * (-1 * BH2_mass * G) / \
        (np.linalg.norm(star_x_vec - BH2_x_vec) ** 3)

    acc_star = acc_star_1 + acc_star_2
    vel_star = star_v_vec

    kwargs['bh1'] = BH1_x_vec
    kwargs['bh2'] = BH2_x_vec

    return np.concatenate((vel_star, acc_star))


def print_default():
    print('\nThese are the default parameters:')
    print('\nDefault X position of Star:\t\t\t\t2.0')
    print('Default Y position of Star:\t\t\t\t0.0')
    print('Default Z position of Star:\t\t\t\t0.0')
    print('Default X component of Velocity of the Star:\t\t0.0')
    print('Default Y component of Velocity of the Star:\t\t1.0')
    print('Default Z component of Velocity of the Star:\t\t0.0')
    print('Default time step: \t\t\t\t\t1.0e-2')
    print('Default maximum run time: \t\t\t\t60')
    print('Default mass ratio: \t\t\t\t\t1.0\n')


def print_help():
    print('\nThis program simulates a neutron star orbiting a binary black hole system')
    print('\npython3 binaryNewton.py')
    print('\n--star\t\t\t\tFlag to begin changing star parameters')
    print('\t--x0,\t-x\t\tSets the x coordinate position of the star')
    print('\t--y0,\t-y\t\tSets the y coordinate position of the star')
    print('\t--z0,\t-z\t\tSets the z coordinate position of the star')
    print('\t--vx0,\t-vx\t\tSets the Vx velocity component of the star')
    print('\t--vy0,\t-vy\t\tSets the Vy velocity component of the star')
    print('\t--vz0,\t-vz\t\tSets the Vz velocity component of the star')
    print('--tstep, -ts\t\t\tSets the time step for the simulation data')
    print('--tmax, -tm\t\t\tSets the maximum run time for the simulation data')
    print('--mratio, -q\t\t\tSets the mass ratio for the binary system')
    print('--sep, -s\t\t\tSets the separation distances of the black holes')
    print('--default, -d\t\t\tShows the default parameters')
    print('--record, -r\t\t\tPrints the initial conditions as a comment\n')


def main(argv):
    # The Initial condition for the star orbiting a black hole
    # The black holes' collective mass is given by the 'mass',
    # Newton's constant by 'G', the star's initial position by
    # (x0, y0, z0), the star's initial velocity by (vx, vy, vz)
    # and q is the mass ratio between the black holes.

    # Initializing default parameters

    # Default mass, G, and q values
    kwargs = {'mass': 1.0, 'G': 1.0, 'q': 1.0}

    x0 = 2.0
    y0 = 0.0
    z0 = 0.0

    vx0 = 0.0
    vy0 = 1.0
    vz0 = 0.0

    # dt is the timestep. The error will be proportional to dt**4
    dt = 1.0e-2

    # start time
    t = 0.0

    # max time
    tmax = 60

    # Black Hole 1's initial position
    BH1x = 1.0
    BH1y = 0.0
    BH1z = 0.0

    # Black Hole 2's initial position
    BH2x = -1.0
    BH2y = 0.0
    BH2z = 0.0

    # Processing command line arguments
    # This will possibly change some of the default values

    i = 0
    record_comment = False

    if len(argv) == 0:
        print('# Running with default settings')

    while i < len(argv):  # while there are unprocessed arguments
        if argv[i] == '--star':
            i += 1
            while i < len(argv):  # while there are unprocessed star
                #print('Star arguments')
                if argv[i] == '--x0' or argv[i] == '-x':
                    i += 1
                    x0 = float(argv[i])
                    #print('X position  changed')
                elif argv[i] == '--y0' or argv[i] == '-y':
                    i += 1
                    y0 = float(argv[i])
                    #print('Y position  changed')
                elif argv[i] == '--z0' or argv[i] == '-z':
                    i += 1
                    z0 = float(argv[i])
                    #print('Z position  changed')
                elif argv[i] == '--vx0' or argv[i] == '-vx':
                    i += 1
                    vx0 = float(argv[i])
                    #print('Velocity x vector  changed')
                elif argv[i] == '--vy0' or argv[i] == '-vy':
                    i += 1
                    vy0 = float(argv[i])
                    #print('Velocity y vector  changed')
                elif argv[i] == '--vz0' or argv[i] == '-vz':
                    i += 1
                    vz0 = float(argv[i])
                    #print('Velocity z vector  changed')
                else:
                    # If the *current* argument is not for a star, counter the *next* increment
                    i -= 1
                    break
                # move to the next argument
                i += 1
        elif argv[i] == '--tstep' or argv[i] == '-ts':
            i += 1
            dt = float(argv[i])
            #print('Time step changed')
        elif argv[i] == '--tmax' or argv[i] == '-tm':
            i += 1
            tmax = float(argv[i])
            #print('Maximum run time changed')
        elif argv[i] == '--mratio' or argv[i] == '-q':
            i += 1
            kwargs['q'] = float(argv[i])
            #print('Mass ratio changed')
        elif argv[i] == '--sep' or argv[i] == '-s':
            i += 1
            BH1x = float(argv[i])/2
            BH2x = -1.0*BH1x
            # Print('Seperation distance of black holes changes)
        elif argv[i] == '--help' or argv[i] == '-h':
            print_help()
            exit(0)
        elif argv[i] == '--default' or argv[i] == '-d':
            print_default()
            exit(0)
        elif argv[i] == '-r' or argv[i] == '--record':
            record_comment = True
        else:
            print('\n "', argv[i], '" is not an option!!')
            print_help()
            exit(1)
        i += 1

    # Calculate initial star position
    initial_position = np.array((x0, y0, z0), dtype=np.float64)

    # Calculate initial star velocity
    initial_velocity = np.array((vx0, vy0, vz0), dtype=np.float64)

    # Concatanate star parameters
    Y = np.concatenate((initial_position, initial_velocity))

    # Puts blach hole 1 position parameters into an array
    initial_bh1_pos = np.array((BH1x, BH1y, BH1z), dtype=np.float64)

    BH1 = initial_bh1_pos

    # Puts black hole 2 position parameters into an array
    initial_bh2_pos = np.array((BH2x, BH2y, BH2z), dtype=np.float64)

    BH2 = initial_bh2_pos

    kwargs['bh1'] = BH1
    kwargs['bh2'] = BH2

    omega, BH_dist = calc_omega(kwargs['mass'], kwargs['G'], BH1, BH2)
    kwargs['omega'] = omega
    kwargs['BH_dist'] = BH_dist

    if record_comment:
        print('# Star Position: x:', x0, ' y:', y0, ' z:', z0)
        print('# Star Velocity Components: vx0: ',
              vx0, ' vy0:', vy0, ' vz0:', vz0)
        print('# Time Step:', dt, '\tRun Time max:', tmax)
        print('')

    """
    list_size = int(tmax / dt) + 1

    r_list = np.zeros(shape=(list_size,))
    phi_list = np.zeros(shape=(list_size,))
    x_list = np.zeros(shape=(list_size,))
    y_list = np.zeros(shape=(list_size,))
    t_list = np.zeros(shape=(list_size,))
    """

    while t < tmax:
        """
        x_pos = Y[0]
        y_pos = Y[1]
        z_pos = Y[2]
        r = np.sqrt(x_pos ** 2 + y_pos ** 2)
        phi = np.arctan2(y_pos, x_pos)
        phi_list[lst_indx] = phi
        r_list[lst_indx] = r
        x_list[lst_indx] = Y[0]
        y_list[lst_indx] = Y[1]
        t_list[lst_indx] = t
        lst_indx += 1
        """

        BH1 = kwargs['bh1']
        BH2 = kwargs['bh2']
        print(t, Y[0], Y[1], Y[2], BH1[0], BH1[1],
              BH1[2], BH2[0], BH2[1], BH2[2])

        # The Runge-Kutta routine returns the new value of Y, t, and a
        # possibly updated value of dt
        t, Y, dt = RK4_Step(t, Y, dt, Keppler_Binary_RHS, **kwargs)


if __name__ == "__main__":
    main(sys.argv[1:])
