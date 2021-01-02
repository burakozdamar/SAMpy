#!/usr/bin/env/python3

import sys
import numpy as np

sys.path.append('..')

ff = sys.argv[1]
gg = "pos_recentered_py.xyz"

from utils.xyz_utils import read_xyz
from collections import Counter
from itertools import cycle 

mass = {'H': 1.00794,'x':0,'C':12.0107,'O':15.9994,'Si':28.0855,'Cl':35.4527,'K':39.0983,'Al':26.981539}
pbc_ = np.array([13.386,13.286,85])


def memoize_mass(ff):
  mass_arr = []
  with open(ff) as f:
    xyz = read_xyz(f)
    atomtypes = xyz.atomtypes
  # silicons or alums coming later from C
  if 'C' in atomtypes:
    mass_arr = np.array([mass[atom] if atom == 'Si' else 0 for atom in atomtypes])
  return mass_arr, atomtypes

mass_arr, atomtypes = memoize_mass(ff)

def translate(arr1,arr2):
  return arr1-arr2

def calcul_CM(xyz):
  coords = xyz.coords
  atomtypes = xyz.atomtypes
  com = np.average(coords, axis=0, weights=mass_arr)
  t = pbc(translate(coords, com))
  return t

def pbc(arr):
  arr = np.where(arr < -pbc_/2, arr+pbc_, arr)
  arr = np.where(arr >  pbc_/2, arr-pbc_, arr)
  return arr 

def write_xyz(fout, coords, title="", atomtypes=("A",)):
    with open(gg,'a') as fout:
      fout.write("%d\n%s\n" % (coords.size / 3, title))
      for x, atomtype in zip(coords.reshape(-1, 3), cycle(atomtypes)):
          fout.write("%s %.8g %.8g %.8g\n" % (atomtype, x[0], x[1], x[2]))

def open_file(ff):

  with open(ff) as f:
    with open(gg,'w') as g:
      for i in range(20):
        try:
          xyz = read_xyz(f)
          print(f"{i}")
          t = calcul_CM(xyz)
          write_xyz(g, t, title="", atomtypes=atomtypes) 
          #write_to_file(g,t)
        except Exception as e:
          print("DONE",e)
          break
  
open_file(ff)

 
