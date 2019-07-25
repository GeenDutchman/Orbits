#!/usr/bin/python3

import numpy as np
from os import walk, getcwd
import sys

print(getcwd())
f = []
for (dirpath, dirnames, filenames) in walk(getcwd()):
    f.extend(filenames)
    break # will only get the files right there

print(f)
sep_set = set([])

Pre_N = input("Enter Pre N ")
Post_N = input("Enter Post N ")
Pre_R = input("Enter Pre R ")
Post_R = input("Enter Post R ")

for filename in f:
    if (Pre_N in filename and Post_N in filename):
        sep_set.add(filename[len(Pre_N):-len(Post_N)])
    elif (Pre_R in filename and Post_R in filename):
        sep_set.add(filename[len(Pre_R):-len(Post_R)])
print(sep_set)

# stop = input("stoppit?q/n ")
# if stop is 'q':
#     exit(0)


sep_list = []
postfix = ""
for bh_sep in sep_set:
    Rfile = Pre_R + bh_sep + Post_R
    Nfile = Pre_N + bh_sep + Post_N
    print(Rfile, Nfile)
    Rdata = np.genfromtxt(Rfile, names=True, usecols=('Precession'))
    Ndata = np.genfromtxt(Nfile, names=True, usecols=('Precession'))

    Rcol = Rdata['Precession']
    Ncol = Ndata['Precession']

    print(Rcol, type(Rcol))
    print(Ncol, type(Ncol))

    Nfloat = isinstance(Ncol, (float, np.ndarray))
    Rfloat = isinstance(Rcol, (float, np.ndarray))

    sep_list.append([float(bh_sep), Ncol if Nfloat else Ncol[0], Rcol if Rfloat else Rcol[0]])

sep_list = sorted(sep_list)
print(sep_list)

header_text = "Separation Newtonian Relativistic"
np.savetxt('result.dat', sep_list, header=header_text)
