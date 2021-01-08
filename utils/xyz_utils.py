import sys
import numpy as np
from collections import namedtuple
from itertools import groupby

# inspired from pele/utils/xyz.py
def read_xyz(f):
    natoms = int(f.readline())

    data = f.readline()[:-1].split()
    step = int(data[2].split(",")[0])

    try:
      time = float(data[5].split(",")[0])
      energy = float(data[8].split(",")[0])
    except:
      time, energy = None, None

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
  print("Reading BOXDATA")
  with open("../warehouse/BOXDATA") as f:
    ll = [l.strip() for l in f.readlines()]
    l = [list(sub) for ele, sub in groupby(ll, key = bool) if ele] 
    k = [item[0] for item in l]
    v = [item[1:] for item in l]
    v = [[i.split(' ') if len(i.split(' '))>1 else i for i in sl] for sl in v]
    d = dict(zip(k,v))
    for k, v in d.items():
      if isinstance(v[0],list):
        d[k]=v[0]
      elif v[0].isdigit():
        d[k]=int(v[0])
      elif 'd' in v[0] and v[0].split('d')[0].lstrip('-').isdigit():
        a = float(v[0].split('d')[0])
        d[k]=a
    return d
