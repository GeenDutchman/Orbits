from __future__ import print_function
import numpy as np
from RK import RK4_Step, RK45_Step
import sys
from PN import center_of_mass_coordinates_to_BH_positions, Omega_of_r, PN_Acceleration

'''
    Helps manage indexing into the Y list by dynamically assigning (key: Y_index) pairs.
    It can also add a list or array of these pairs given two lists of identical length.
'''


def addY(values, keys, Y, key_dict):
    if isinstance(values, (list, tuple, np.ndarray)) and isinstance(keys, (list, tuple, np.ndarray)):
        if len(values) != len(keys):
            raise IndexError(
                'The length of values must match the length of keys!!')
        for index in range(len(values)):
            Y_len = len(Y)
            Y = np.append(Y, values[index])
            key_dict[keys[index]] = Y_len
    else:
        Y_len = len(Y)
        Y = np.append(Y, values)
        key_dict[keys] = Y_len
    return Y, key_dict


'''
    This method is responsible for calculating the time derivatives of each thing in y0.
    Additional arguments are passed to subordinate methods through kwargs, most of which define
    default values if necessary.
'''


def Keppler_Binary_RHS(t, y0, **kwargs):
    if 'Y_dict' in kwargs:
        Y_dict = kwargs['Y_dict']
    else:
        print('Must have the dictionary!!')
        exit(2)

    star_x_vec = np.array([y0[Y_dict['star_x']],
                           y0[Y_dict['star_y']], y0[Y_dict['star_z']]], dtype=np.float64)  # star position
    star_v_vec = np.array([y0[Y_dict['star_vx']],
                           y0[Y_dict['star_vy']], y0[Y_dict['star_vz']]], dtype=np.float64)  # star velocity
    # star angle position
    #star_phi_vec = y0[Y_dict['star_angle']]

    # star position vector norm, such that
    # star_r = (star_x**2 + star_y**2 + star_z**2)**1/2
    # or the distance from the origin
    # star_r = np.linalg.norm(star_x_vec)

    '''
        This part used to have a bunch of "if 'this' in kwargs, set it from kwargs. Else, here is some default value.
        Most of those are found in subordinate methods now.
    '''

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

    # corrections: star_acc, bh_r_dot, bh_Omega_dot
    star_v_dot, bh_r_dot, bh_Omega_dot = PN_Acceleration(
        star_x_vec, star_v_vec, y0[Y_dict['bh_r']], y0[Y_dict['bh_psi']], y0[Y_dict['bh_Omega']], **kwargs)

    # calculate the current position, but does not do the z coord??
    BH1_x_vec[0:3], BH2_x_vec[0:3] = center_of_mass_coordinates_to_BH_positions(
        y0[Y_dict['bh_r']], y0[Y_dict['bh_psi']], **kwargs)

    # One way to calculate mass of the Black Holes
    # BH1_mass = combined_BH_mass / (BH_ratio + 1)
    # BH2_mass = (-1 * combined_BH_mass * BH_ratio) / (BH_ratio + 1)

    # star angle
    phi_dot = np.linalg.norm(
        np.cross(star_x_vec, star_v_vec)) / (np.linalg.norm(star_x_vec) ** 2)

    kwargs['bh1'] = BH1_x_vec
    kwargs['bh2'] = BH2_x_vec

    # holds all the changes
    deltas = np.array([None] * len(y0))

    # put the dot where the original goes
    deltas[Y_dict['star_x']] = star_v_vec[0]
    deltas[Y_dict['star_y']] = star_v_vec[1]
    deltas[Y_dict['star_z']] = star_v_vec[2]

    deltas[Y_dict['star_vx']] = star_v_dot[0]
    deltas[Y_dict['star_vy']] = star_v_dot[1]
    deltas[Y_dict['star_vz']] = star_v_dot[2]

    deltas[Y_dict['star_angle']] = phi_dot

    deltas[Y_dict['bh_r']] = bh_r_dot
    deltas[Y_dict['bh_psi']] = y0[Y_dict['bh_Omega']]  # bh_psi_dot is bh_Omega
    deltas[Y_dict['bh_Omega']] = bh_Omega_dot

    return deltas


'''
    Simple enough, this will print the default values for the various arguments.
'''


def print_default():
    print('\nThese are the default parameters:')
    print('\nDefault X position of Star:\t\t\t\t5000.0')
    print('Default Y position of Star:\t\t\t\t0.0')
    print('Default Z position of Star:\t\t\t\t0.0')
    print('Default X component of Velocity of the Star:\t\t0.0')
    print('Default Y component of Velocity of the Star:\t\t0.014')
    print('Default Z component of Velocity of the Star:\t\t0.0')
    print('Default time step: \t\t\t\t\t1.0e-2')
    print('Default maximum run time: \t\t\t\t1.0e12')
    print('Default black hole separation: \t\t\t\t100')
    print('Default mass ratio: \t\t\t\t\t1.0')
    print('Default maximum orbits: \t\t\tinfinite\n')


'''
    This prints the help
'''


def print_help():
    print('\nThis program simulates a neutron star orbiting a binary black hole system')
    print('\npython3 relativisticBinary.py')
    print('\n--star\t\t\t\tFlag to begin changing star parameters')
    print('\t--x0,\t-x\t\tSets the x coordinate position of the star')
    print('\t--y0,\t-y\t\tSets the y coordinate position of the star')
    print('\t--z0,\t-z\t\tSets the z coordinate position of the star')
    print('\t--vx0,\t-vx\t\tSets the Vx velocity component of the star')
    print('\t--vy0,\t-vy\t\tSets the Vy velocity component of the star')
    print('\t--vz0,\t-vz\t\tSets the Vz velocity component of the star')
    print('--tstep, -ts\t\t\tSets the time step for the simulation data')
    print('--tmax, -tm\t\t\tSets the maximum run time for the simulation data')
    print('--omax, -om\t\t\tSets the maximum number of orbits for the simulation data')
    print('--mratio, -q\t\t\tSets the mass ratio for the binary system')
    print('--sep, -s\t\t\tSets the separation distances of the black holes')
    print('\t--rx,\t-x\t\tSets the x componet of the black hole separation')
    print('\t--ry,\t-y\t\tSets the y componet of the black hole separation')
    print('--default, -d\t\t\tShows the default parameters')
    print('--extended, -e\t\t\tWill print extra data to the file')
    print('--rk45, -45\t\t\tSets to auto adjust the time-step dynamically\n')


'''
    This is more for clerical information when plotting.

    By keeping track of the star's orbit, all bodies can be kept in view.
'''


def update_min_max(min_max_list, Y, index):
    # if new coordinate is less than stored minimum, update
    if Y[index] < min_max_list[0]:
        min_max_list[0] = np.floor(Y[index])
    elif Y[index] > min_max_list[1]:  # if new coord is more than stored maximum
        min_max_list[1] = np.ceil(Y[index])


'''
    Main set up default values, parses command line arguments, and then if parameters are
    OK, will run the simulation until A) the simulation time is elapsed, B) the orbital maximum is reached, 
    or C) the Black Holes are at or within 10 separations of eachother.
'''


def main(argv):
    # The Initial condition for the star orbiting a black hole
    # The black holes' collective mass is given by the 'mass',
    # Newton's constant by 'G', the star's initial position by
    # (x0, y0, z0), the star's initial velocity by (vx, vy, vz)
    # and q is the mass ratio between the black holes.

    # Initializing default parameters

    # Default mass, G, and q values
    kwargs = {'mass': 1.0, 'G': 1.0, 'massratio': 1.0}

    # Default star's position
    x0 = 5000
    y0 = 0.0
    z0 = 0.0
    # Default star's velocity
    vx0 = 0.0
    vy0 = 0.014
    vz0 = 0.0

    # dt is the timestep. The error will be proportional to dt**4
    dt = 1.0e-2

    # start time
    t = 0.0

    # max time
    tmax = 1.0e0300

    # max orbits, if this is '-1' the orbital maximum will be ignored
    MAX_ORBITS = -1

    # Separation of Black Holes' initial position
    r_x_hat = 100.0
    r_y_hat = 0.0
    r_z_hat = 0.0

    # Processing command line arguments
    # This will possibly change some of the default values

    i = 0
    use_RK_45 = False
    extended = 0  # controls how much information is printed

    kwargs['tol'] = 1e-8  # tolerance for RK_45

    # 'R' relativistic '1e2' separations 'x50' star is 50 times the separation 'T1e-8' tolerance
    file_name = "R1e2x50T1e-8.dat"

    # for better options menu https://docs.python.org/3/library/argparse.html#sub-commands
    # But we just kept it simple

    while i < len(argv):  # while there are unprocessed arguments
        if argv[i] == '--star':
            i += 1
            while i < len(argv):  # while there are unprocessed star arguments
                # print('Star arguments')
                if argv[i] == '--x0' or argv[i] == '-x':
                    i += 1
                    x0 = float(argv[i])
                    # automatically calculated to be less than circular velocity
                    vy0 = 1/(np.sqrt(x0)) * 0.99
                elif argv[i] == '--y0' or argv[i] == '-y':
                    i += 1
                    y0 = float(argv[i])
                elif argv[i] == '--z0' or argv[i] == '-z':
                    i += 1
                    z0 = float(argv[i])
                elif argv[i] == '--vx0' or argv[i] == '-vx':
                    i += 1
                    vx0 = float(argv[i])
                elif argv[i] == '--vy0' or argv[i] == '-vy':
                    i += 1
                    vy0 = float(argv[i])
                elif argv[i] == '--vz0' or argv[i] == '-vz':
                    i += 1
                    vz0 = float(argv[i])
                else:
                    # If the *current* argument is not for a star, counter the *next* increment
                    i -= 1
                    break
                # move to the next argument
                i += 1
        elif argv[i] == '--sep':
            i += 1
            while i < len(argv):  # while there are unprocessed separation arguments
                # print('Star arguments')
                if argv[i] == '--rx' or argv[i] == '-x':
                    i += 1
                    r_x_hat = float(argv[i])
                elif argv[i] == '--ry' or argv[i] == '-y':
                    i += 1
                    r_y_hat = float(argv[i])
                else:
                    # If the *current* argument is not for a separation, counter the *next* increment
                    i -= 1
                    break
                # move to the next argument
                i += 1
        elif argv[i] == '--tstep' or argv[i] == '-ts':
            i += 1
            dt = float(argv[i])
            # print('Time step changed')
        elif argv[i] == '--tmax' or argv[i] == '-tm':
            i += 1
            tmax = float(argv[i])
            # print('Maximum run time changed')
        elif argv[i] == '--omax' or argv[i] == '-om':
            i += 1
            MAX_ORBITS = float(argv[i])
            # print('Maximum number of orbits changed')
        elif argv[i] == '--mratio' or argv[i] == '-q':
            i += 1
            kwargs['massratio'] = float(argv[i])
            # print('Mass ratio changed')
        elif argv[i] == '--help' or argv[i] == '-h':
            print_help()
            exit(0)
        elif argv[i] == '--default' or argv[i] == '-d':
            print_default()
            exit(0)
        elif argv[i] == '--rk45' or argv[i] == '-45':
            use_RK_45 = True
        elif argv[i] == '--file' or argv[i] == '-f':
            i += 1
            file_name = argv[i]
        elif argv[i] == '--tol' or argv[i] == '-to':
            i += 1
            kwargs['tol'] = float(argv[i])
        elif argv[i] == '--extended' or argv[i] == '-e':
            i += 1
            extended = int(argv[i])
        else:
            print('\n "', argv[i], '" is not an option!!')
            print_help()
            exit(1)
        i += 1

    # Calculate initial star position
    initial_position = np.array((x0, y0, z0), dtype=np.float64)

    # Calculate initial star velocity
    initial_velocity = np.array((vx0, vy0, vz0), dtype=np.float64)

    # Calculate initial star angle(phi)
    # TODO: need to check if this is the right initialization
    initial_phi = np.array(np.arctan2(y0, x0), dtype=np.float64)

    # Concatenate star parameters
    Y = []  # contains the things to be integrated
    Y_dict = {}  # uses a key:index pairing to keep track of where what values are in Y
    Y, Y_dict = addY(initial_position, [
                     'star_x', 'star_y', 'star_z'], Y, Y_dict)
    Y, Y_dict = addY(initial_velocity, [
                     'star_vx', 'star_vy', 'star_vz'], Y, Y_dict)
    Y, Y_dict = addY(initial_phi, 'star_angle', Y, Y_dict)

    # Puts separation parameters into an array
    initial_separation = np.array(
        (r_x_hat, r_y_hat, r_z_hat), dtype=np.float64)

    r_vec = initial_separation

    # Puts black hole 1 position parameters into an array
    initial_bh1_pos = (r_vec/(1 + kwargs['massratio']))

    BH1 = initial_bh1_pos

    # Puts black hole 2 position parameters into an array
    initial_bh2_pos = (
        (-1 * kwargs['massratio'])/(1 + kwargs['massratio']))*r_vec

    BH2 = initial_bh2_pos

    kwargs['bh1'] = BH1
    kwargs['bh2'] = BH2

    # omega, r_vec = calc_omega(kwargs['mass'], kwargs['G'], BH1, BH2)
    # kwargs['omega'] = omega
    kwargs['BH_dist'] = r_vec

    try:
        f = open(file_name, "x")    # Open a file

        bh_r = np.linalg.norm(r_vec)
        if bh_r <= 10:
            print(
                "# The magnitude of the separation must be larger than 10 separations!!", file=f)
            exit(2)
        Omega = Omega_of_r(bh_r, **kwargs)
        psi = 0

        Y, Y_dict = addY([bh_r, psi, Omega], [
                         'bh_r', 'bh_psi', 'bh_Omega'], Y, Y_dict)

        if extended == 0:
            print('#', 'time', 'star_x', 'star_angle', 'star_r', 'bh_r', file=f)
        elif extended >= 1:
            # use this option if you want to use any of our '.plt' scripts
            print('#', 'time', 'star_x', 'star_y', 'star_z', 'bh1_x', 'bh1_y', 'bh1_z',
                  'bh2_x', 'bh2_y', 'bh2_z', 'star_r', 'star_angle', 'bh_r', 'star_r_dot', end=' ', file=f)
            if extended >= 2:
                # more data
                print('star_vx', 'star_vy', 'star_vz',
                      'bh_psi', 'bh_Omega', end=' ', file=f)
            print('', file=f)

        star_x_min_max = [Y[Y_dict['star_x']], Y[Y_dict['star_x']]]
        star_y_min_max = [Y[Y_dict['star_y']], Y[Y_dict['star_y']]]
        star_z_min_max = [Y[Y_dict['star_z']], Y[Y_dict['star_z']]]

        kwargs['Y_dict'] = Y_dict

        kwargs['use_newtonian_acceleration_in_PN'] = True
        kwargs['use_newtonian_velocities_in_PN'] = False
        # Note that this will override use_newtonian_velocities_in_PN and use_newtonian_acceleration_in_PN if set to False
        kwargs['use_post_newtonian_corrections_in_PN'] = True

        # while time less than max and bh_separation more than 10
        while t < tmax and Y[Y_dict['bh_r']] > 10:

            # star position
            star_pos = [Y[Y_dict['star_x']],
                        Y[Y_dict['star_y']], Y[Y_dict['star_z']]]
            # star velocity
            star_vel = [Y[Y_dict['star_vx']],
                        Y[Y_dict['star_vy']], Y[Y_dict['star_vz']]]

            # Star's distance from the origin
            star_r = np.linalg.norm(star_pos)

            # Dot product are velocity and position of star divided by
            # norm of star position vector
            star_r_dot = (np.dot(star_vel, star_pos))/star_r

            BH1 = kwargs['bh1']
            BH2 = kwargs['bh2']

            # Prints out Time, Star x position, Star y position, Star z Position
            # Black hole 1 x position, Black hole 1 y position, Black hole 1 z position
            # Black hole 2 x position, Black hole 2 y position, Black hole 2 z position
            # Stars distance from origin, Star theta angle relative to origin
            # Black holes' separation distance from each other
            # r(star position) as it changes with respect to time
            # Star x velocity, Star y velocity, Star z velocity

            if extended == 0:
                print(t, Y[Y_dict['star_x']], Y[Y_dict['star_angle']],
                      star_r, Y[Y_dict['bh_r']], file=f)
            elif extended >= 1:
                print(t, Y[Y_dict['star_x']], Y[Y_dict['star_y']], Y[Y_dict['star_z']], BH1[0], BH1[1],
                      BH1[2], BH2[0], BH2[1], BH2[2], star_r, Y[Y_dict['star_angle']], Y[Y_dict['bh_r']], star_r_dot, end=' ', file=f)
                if extended >= 2:
                    print(Y[Y_dict['star_vx']], Y[Y_dict['star_vy']], Y[Y_dict['star_vz']],
                          Y[Y_dict['bh_psi']], Y[Y_dict['bh_Omega']], end=' ', file=f)
                print(file=f)

            # The Runge-Kutta routine returns the new value of Y, t, and a
            # possibly updated value of dt
            if use_RK_45:
                # fine-tunes the dt
                t, Y, dt = RK45_Step(t, Y, dt, Keppler_Binary_RHS, **kwargs)
            else:
                # does not change the dt
                t, Y, dt = RK4_Step(t, Y, dt, Keppler_Binary_RHS, **kwargs)

            # Keeps track of the min and max of the star for plotting purposes
            update_min_max(star_x_min_max, Y, Y_dict['star_x'])
            update_min_max(star_y_min_max, Y, Y_dict['star_y'])
            update_min_max(star_z_min_max, Y, Y_dict['star_z'])

            # only do MAX_ORBITS of...well...orbits
            if MAX_ORBITS > 0 and Y[Y_dict['star_angle']] / (2 * np.pi) >= MAX_ORBITS:
                print('# Maximum orbits: ', MAX_ORBITS, 'reached!!', file=f)
                break

        print('# The star does', Y[Y_dict['star_angle']
                                   ] / (2 * np.pi), 'orbits.', file=f)
        print("# Xmin\tXmax\tYmin\tYmax\tZmin\tZmax", file=f)
        print("#", star_x_min_max[0], "\t", star_x_min_max[1], "\t", star_y_min_max[0],
              "\t", star_y_min_max[1], "\t", star_z_min_max[0], "\t", star_z_min_max[1], file=f)

        f.close()
    except FileExistsError:
        print('A file already exists with that data!!')


if __name__ == "__main__":
    main(sys.argv[1:])
