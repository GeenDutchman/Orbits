#!/usr/bin/python3

import numpy as np
import subprocess as sbp
import sys, os

time_out = 60 * 10 # 10 minute timeout
exponent_min = -5

def relativistic_func(bh_sep, scale, stats):
    # STARTS a separate process runing grScript, and then counts the error if any when the process is done

    exponet = -9

    star_dist = bh_sep * scale  # double check this scaling

    stats['rel_tests'] += 1

    while True:
        try:
            Rfile = "./relativistic/R" + \
                "{:.1e}".format(bh_sep) + "x" + str(scale) + \
                'T' + str(exponet) + ".analyzer.dat"
            # This will run as a separate process, and will timeout beyond time_out
            Rresult = sbp.run(
                ["./grScript.sh", "--star", "-x",  str(star_dist), "-45", "--omax", "10", "--tmax", "1e20", "-ts", "1", "--sep", "-x", str(bh_sep), "-f", Rfile, "--tol", str(1.0 * 10.0 ** exponet)], stdout=sbp.PIPE, stderr=sbp.PIPE, universal_newlines=True, timeout=time_out)
        except TimeoutError as e:
            print('This run went for too long...')
            print(e)
            stats['rel_timeouts'] += 1
            exponet += 1
            if os.path.isfile(Rfile):
                print('Removing file.', Rfile)
                os.remove(Rfile)
            if exponet <= exponent_min:
                print('Trying again with exponet', exponet)
            else:
                print("Cannot go lower than", exponent_min, "as exponet.")
                print()
                return
        else:
            if Rresult.returncode != 0:
                print('There was an error!!', file=sys.stderr)
                stats['rel_err_count'] += 1
                print(Rresult, file=sys.stderr, flush=True)
            else:
                print('Ran successfully...check the output.')
                print(Rresult)
            print()
            return

def newtonian_func(bh_sep, scale, stats):
    # STARTS a separate process running grScript, and then counts the error if any when the process is done

    exponet = -9

    star_dist = bh_sep * scale  # double check this scaling
    star_y_vel = 1 / np.sqrt(star_dist) # for a circular orbit

    stats['newt_tests'] += 1

    while True:
        try:
            Nfile = "./newtonian/N" + \
                "{:.1e}".format(bh_sep) + "x" + str(scale) + \
                'T' + str(exponet) + ".analyzer.dat"
            # runs as a separate process and will timeout after time_out
            Nresult = sbp.run(
                ["./newtonianScript.sh", "--star", "-x", str(star_dist), "-vy", str(star_y_vel), "-45", "--omax", "10", "--tmax", "1e20", "-ts", "1", "--sep", str(bh_sep), "-f", Nfile, "--tol", str(1.0 * 10.0 ** exponet)], stdout=sbp.PIPE, stderr=sbp.PIPE, universal_newlines=True, timeout=time_out)
        except TimeoutError as e:
            print('This run went for too long...')
            print(e)
            stats['newt_timeouts'] += 1
            exponet += 1
            if os.path.isfile(Nfile):
                print('Removing file.', Nfile)
                os.remove(Nfile)
            if exponet <= exponent_min:
                print('Trying again with exponet', exponet)
            else:
                print("Cannot go lower than", exponent_min, "as exponet.")
                print()
                return
        else:
            if Nresult.returncode != 0:
                print('There was an error!!', file=sys.stderr)
                stats['newt_err_count'] += 1
                print(Nresult, file=sys.stderr, flush=True)                
            else:
                print('Ran successfully...check the output')
                print(Nresult)
            print()
            return
        

# this is meant to run a series of values
# as a process in the background
# also make sure to not run any plates or the movie

stats = { "rel_tests": 0, "newt_tests": 0, "rel_err_count": 0, "newt_err_count": 0, "sep_tested": {}, "newt_timeouts": 0, "rel_timeouts": 0}

start_bh_sep_dist = 100
max_sep_dist = 1000
sep_increments = 10 # multiplies, does not add

scale_start = 10
max_scale = 10000
scale_increments = 10 # multiplies, does not add

bh_sep = start_bh_sep_dist
scale = scale_start

# for paralell ideas: https://stackoverflow.com/a/23616229

# increments in terms of a tenth of the base, up to 2 * the base
while bh_sep < max_sep_dist + 1:
    scale = scale_start
    stats['sep_tested'][str(bh_sep)] = []
    while scale < max_scale + 1:
        stats['sep_tested'][str(bh_sep)].append(scale)
        relativistic_func(bh_sep, scale, stats)
        newtonian_func(bh_sep, scale, stats)

        # now increment scale
        scale *= scale_increments

    # now increment separation
    bh_sep *= sep_increments

# Alert the user that the code has finished
msg_string = 'The test is complete.\n' + str(stats)

result = sbp.run(["which", "mail"], stdout=sbp.DEVNULL)
if result.returncode == 0:
    msg = sbp.Popen(['printf', msg_string], stdout=sbp.PIPE)
    result = sbp.run(["mail", "-s", "Test run", "GeenDutchman@mail.rit.edu", "larreu@rit.edu"], stdin=msg.stdout)
    msg.wait()
result = sbp.run(["which", "espeak"], stdout=sbp.DEVNULL)
if result.returncode == 0:
    result = sbp.run(["espeak", msg_string])
print(msg_string)
# if neither goes off....oh well
