#!/usr/bin/python3

import numpy as np
import subprocess as sbp

# this is meant to run a series of values
# as a process in the background
# also make sure to not run any plates or the movie

star_dist = 200
star_y_vel = 0.05
base_bh_sep = 200

# increments in terms of a tenth of the base, up to 2 * the base
for bh_sep in range(base_bh_sep, base_bh_sep * 5, int(base_bh_sep * 0.1)):
    star_dist = bh_sep * 50 # double check this scaling
    star_y_vel = 1 / np.sqrt(star_dist) # for a circular orbit
    result = sbp.run(
        ["./grScript.sh", "--star", "-x",  str(star_dist), "-vy", str(star_y_vel), "-45", "--omax", "100", "--tmax", "1e12", "--sep", "-x", str(bh_sep), "-f", "./relativistic/R" + str(bh_sep) + ".dat"])

# Alert the user that the code has finished
result = sbp.run(["which", "mail"], stdout=sbp.DEVNULL)
if result.returncode == 0:
    msg = sbp.Popen(['echo', 'The test is complete.'], stdout=sbp.PIPE)
    result = sbp.run(["mail", "-s", "Test run", "GeenDutchman@mail.rit.edu", "larreu@rit.edu"], stdin=msg.stdout)
    msg.wait()
result = sbp.run(["which", "espeak"], stdout=sbp.DEVNULL)
if result.returncode == 0:
    result = sbp.run(["espeak", "This run is complete."])
# print(result)
# if neither goes off....oh well
