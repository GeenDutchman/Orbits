#!/usr/bin/python3

import numpy as np
import subprocess as sbp

for i in range(100, 102, 1):
    # this is meant to run a series of values
    # as a process in the background
    # also make sure to not run any plates or the movie
    result = sbp.run(
        ["./grScript.sh", "--star", "-x",  "200", "-vy", "0.05" ,"-45", "--tmax", "9e6", "--sep", "-x", (i).__str__()]) #, stdout=sbp.DEVNULL, stderr=sbp.DEVNULL)

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
