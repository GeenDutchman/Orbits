#!/usr/bin/python3

import numpy as np
import subprocess as sbp
import sys

time_out = 60 * 15 # 15 minute timeout
exponent_min = -5

def relativistic_func(bh_sep, star_x, star_vel, stats):
    # STARTS a separate process runing grScript, and then counts the error if any when the process is done

    exponet = -9

    Rfile = "./relativistic/R" + str(bh_sep) + ".dat"
    stats['rel_tests'] += 1

    while True:
        try:
            Rresult = sbp.run(
                ["./grScript.sh", "--star", "-x",  str(star_dist), "-vy", str(star_y_vel), "-45", "--omax", "100", "--tmax", "1e12", "-tm", "1", "--sep", "-x", str(bh_sep), "-f", Rfile, "--tol", str(1.0 * 10.0 ** exponet)], stdout=sbp.PIPE, stderr=sbp.PIPE, universal_newlines=True)
        except TimeoutError as e:
            print('This run went for too long...')
            print(e)
            stats['rel_timeouts'] += 1
            exponet += 1
            if exponet <= exponent_min:
                print('Trying again with exponet', exponet)
            else:
                print("Cannot go lower than", exponent_min, "as exponet.")
                return
        else:
            if Rresult.returncode != 0:
                print('There was an error!!', file=sys.stderr)
                stats['rel_err_count'] += 1
                print(Rresult, file=sys.stderr, flush=True)
            else:
                print('Ran successfully...check the output.')
            return

def newtonian_func(bh_sep, star_x, star_vel, stats):
    # STARTS a separate process running grScript, and then counts the error if any when the process is done

    exponet = -9


    Nfile = "./newtonian/N" + str(bh_sep) + ".dat"
    stats['newt_tests'] += 1

    while True:
        try:
            Nresult = sbp.run(
                ["./newtonianScript.sh", "--star", "-x", str(star_dist), "-vy", str(star_y_vel), "-45", "--omax", "100", "--tmax", "1e12", "-tm", "1", "--sep", str(bh_sep), "-f", Nfile, "--tol", str(1.0 * 10.0 ** exponet)], stdout=sbp.PIPE, stderr=sbp.PIPE, universal_newlines=True)
        except TimeoutError as e:
            print('This run went for too long...')
            print(e)
            stats['newt_timeouts'] += 1
            exponet += 1
            if exponet <= exponent_min:
                print('Trying again with exponet', exponet)
            else:
                print("Cannot go lower than", exponent_min, "as exponet.")
                return
        else:
            if Nresult.returncode != 0:
                print('There was an error!!', file=sys.stderr)
                stats['newt_err_count'] += 1
                print(Nresult, file=sys.stderr, flush=True)                
            else:
                print('Ran successfully...check the output')
            return
        

# this is meant to run a series of values
# as a process in the background
# also make sure to not run any plates or the movie

stats = { "rel_tests": 0, "newt_tests": 0, "rel_err_count": 0, "newt_err_count": 0, "sep_tested": [], "newt_timeouts": 0, "rel_timeouts": 0}

start_bh_sep_dist = 100
max_sep_dist = 100#1000
sep_increments = 10 # multiplies, does not add

scale_start = 10
max_scale = 100#10000
scale_increments = 10 # multiplies, does not add

bh_sep = start_bh_sep_dist
scale = scale_start

# for paralell ideas: https://stackoverflow.com/a/23616229

# increments in terms of a tenth of the base, up to 2 * the base
while bh_sep < max_sep_dist + 1:
    scale = scale_start
    while scale < max_scale + 1:

        star_dist = bh_sep * scale # double check this scaling
        star_y_vel = 1 / np.sqrt(star_dist) # for a circular orbit
        stats['sep_tested'].append(bh_sep)
        
        relativistic_func(bh_sep, star_dist, star_y_vel, stats)
        newtonian_func(bh_sep, star_dist, star_y_vel, stats)

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
# print(result)
# if neither goes off....oh well
