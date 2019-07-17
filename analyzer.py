#!/usr/bin/python3

import numpy as np
import subprocess as sbp
import sys

# this is meant to run a series of values
# as a process in the background
# also make sure to not run any plates or the movie

rel_tests = 0
newt_tests = 0
rel_err_count = 0
newt_err_count = 0
separations_tested = []

start_bh_sep_dist = 200
max_sep_dist = 1000
increments = 100

# for paralell ideas: https://stackoverflow.com/a/23616229

# increments in terms of a tenth of the base, up to 2 * the base
for bh_sep in range(start_bh_sep_dist, max_sep_dist + 1, increments):
    star_dist = bh_sep * 50 # double check this scaling
    star_y_vel = 1 / np.sqrt(star_dist) # for a circular orbit
    separations_tested.append(bh_sep)
    Rresult = sbp.run(
       ["./grScript.sh", "--star", "-x",  str(star_dist), "-vy", str(star_y_vel), "-45", "--omax", "100", "--tmax", "1e12", "--sep", "-x", str(bh_sep), "-f", "./relativistic/R" + str(bh_sep) + ".dat"])
    rel_tests += 1
    if Rresult.returncode != 0:
        print('There was an error!!', file=sys.stderr)
        rel_err_count += 1
        print(Rresult, file=sys.stderr, flush=True)
    Nresult = sbp.run(
        ["./newtonianScript.sh", "--star", "-x", str(star_dist), "-vy", str(star_y_vel), "-45", "--omax", "100", "--tmax", "1e12", "--sep", str(bh_sep), "-f", "./newtonian/N" + str(bh_sep) + ".dat"])
    newt_tests += 1
    if Nresult.returncode != 0:
        print('There was an error!!', file=sys.stderr)
        newt_err_count += 1
        print(Nresult, file=sys.stderr, flush=True)

# Alert the user that the code has finished
msg_string = 'The test is complete.\nThere were ' + str(rel_err_count) + ' relativistic errors out of ' + str(rel_tests) + ' tests and ' + str(newt_err_count) + ' newtonian errors out of ' + str(newt_tests) + ' tests.\nThis test used the separations:\n' + str(separations_tested) + '\n'

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
