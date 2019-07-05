#!/usr/bin/python3

import numpy as np
import subprocess as sbp

for i in range(100, 103, 1):
    # this is meant to run a series of values
    # as a process in the background
    # also make sure to not run any plates or the movie
    result = sbp.run(
        ["./grScript.sh", "--star", "-x",  "200", "-vy", "0.05" ,"-45", "--tmax", "9e6", "--sep -x", (i).__str__()], stdout=sbp.DEVNULL, stderr=sbp.DEVNULL)

# Alert the user that the code has finished
result = sbp.run(["which mail"], shell=True)
if result.returncode == 0:
    result = sbp.run(["mail -s 'Test run' GeenDutchman@mail.rit.edu larreu@rit.edu <<< 'The test is complete'"], shell=True)
result = sbp.run(["which espeak"], shell=True, stdout=sbp.DEVNULL)
if result.returncode == 0:
    result = sbp.run(["espeak", "'This run is complete.'"], shell=True)
print(result)
# if neither goes off....oh well
