#!/usr/bin/python3

import numpy as np
import subprocess


for i in range(0, 42, 2):
    # this is meant to run a series of values
    # as a process in the background
    # also make sure to not run any plates or the movie
    return_code = subprocess.run(
        ["./bashScript.sh", "--star", "-x",  "100", "-vy", "0.0" ,"-vz", "-0.05", "-45", "--tmax", "25000", "--sep", (i*0.1).__str__()], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

