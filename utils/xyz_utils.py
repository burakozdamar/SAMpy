import sys
import numpy as np
from collections import namedtuple
from itertools import groupby

# inspired from pele/utils/xyz.py
def read_xyz(f):
    natoms = int(f.readline())

    data = f.readline()[:-1].split()
    step = int(data[2].split(",")[0])
    time = float(data[5].split(",")[0])
    energy = float(data[8].split(",")[0])

    data = np.array([step, time, energy])
    coords= np.zeros([natoms, 3], dtype="float64")
    atomtypes = []

    for x in coords:
        line = f.readline().split()
        atomtypes.append(line[0])
        x[:] = list(map(float, line[1:4]))

    return namedtuple("XYZFile", ["coords", "data", "atomtypes"])(
        coords, data, atomtypes
    )

def read_boxdata():
  with open("../warehouse/BOXDATA") as f:
    ll = [l.strip() for l in f.readlines()]
    res = [list(sub) for ele, sub in groupby(ll, key = bool) if ele] 
    k = [item[0] for item in res]
    v = [item[1:] for item in res]
    v = [[i.split(' ') if len(i.split(' '))>1 else i for i in sl] for sl in v]
    return dict(zip(k,v))
