import sys
import numpy as np
from collections import namedtuple

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


def hhread_xyz(fin):
    natoms = int(fin.readline())
    step = int(fin.readline()[:-1].split()[-1])
    coords = np.zeros([natoms, 3], dtype="float64")
    atomtypes = []
    for x in coords:
        line = fin.readline().split()
        atomtypes.append(line[0])
        x[:] = list(map(float, line[1:4]))

    return namedtuple("XYZFile", ["coords", "step", "atomtypes"])(
        coords, step, atomtypes
    )

