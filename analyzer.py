#!/usr/bin/python3

import numpy as np
import subprocess


for i in range(0, 4, 2):
    return_code = subprocess.run(
        ["./bashScript.sh", "--star", "-x",  "100", "-vy", "0.05", "-45", "--tmax", "25000", "--sep", (i*0.1).__str__()], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

