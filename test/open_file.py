#!/usr/bin/env/python3

import re
import sys
import numpy as np

sys.path.append('..')

ff = sys.argv[1]
gg = "py_pos_recentered.xyz"
hh = "py_pos_rebuilt.xyz"

from utils.xyz_utils import read_xyz
from collections import Counter, namedtuple
from itertools import cycle, count, groupby

mass = {'H': 1.00794,'x':0,'C':12.0107,'O':15.9994,'Si':28.0855,'Cl':35.4527,'K':39.0983,'Al':26.981539}
pbc_ = np.array([13.386,13.286,85])

def memoize_mass(xyz):
  print('here')
  print(xyz.atomtypes)
  mass_arr = np.array([mass[atom] for atom in xyz.atomtypes])  
  return mass_arr 
 # if kind=='water':
 #   mass_arr = np.array([mass[atom] for atom in atomtypes])
 # else:
 #   if 'C' in atomtypes:
 #     mass_arr = np.array([mass[atom] if atom == 'Si' else 0 for atom in atomtypes])
 #   else:
 #     mass_arr = np.array([mass[atom] for atom in atomtypes])
 # return mass_arr, atomtypes

#mass_arr, atomtypes = memoize_mass(ff)

def translate(arr1,arr2):
  return arr1-arr2

def calcul_CM(xyz, kind='Si'):
  coords = xyz.coords
  
  if kind=='water':
    mass_arr = np.array([mass[atom] for atom in xyz.atomtypes])
    print(mass_arr)
  else:
    mass_arr = np.array([mass[atom] if atom == 'Si' else 0 for atom in xyz.atomtypes])  

  com = np.average(coords, axis=0, weights=mass_arr[:coords.shape[0]])

  t = pbc(translate(coords, com))
  return t

def water_molecules(xyz):
  data = xyz.data
  step = data[0]
  atomtypes = xyz.atomtypes
  coords = xyz.coords
  data_dct = read_boxdata()
  val = data_dct['$NO']
  coords = coords[:int(val[0])*3]  
  atomtypes = atomtypes[:int(val[0])*3]  
  return namedtuple("waters", ["coords", "step", "atomtypes"])(
        coords, step, atomtypes 
    )

#def read_xyz(fin):
#    natoms = int(fin.readline())
#    step = int(fin.readline()[:-1].split()[-1])
#    coords = np.zeros([natoms, 3], dtype="float64")
#    atomtypes = []
#    for x in coords:
#        line = fin.readline().split()
#        atomtypes.append(line[0])
#        x[:] = list(map(float, line[1:4]))
#
#    return namedtuple("XYZFile", ["coords", "step", "atomtypes"])(
#        coords, step, atomtypes
#    )


def pbc(arr):
  arr = np.where(arr < -pbc_/2, arr+pbc_, arr)
  arr = np.where(arr >  pbc_/2, arr-pbc_, arr)
  return arr 

def write_xyz(fout, coords, title="", atomtypes=("A",)):
  fout.write("%d\n%s\n" % (coords.size / 3, title))
  for x, atomtype in zip(coords.reshape(-1, 3), cycle(atomtypes)):
    #fout.write("%s %.8g %.8g %.8g\n" % (atomtype, x[0], x[1], x[2]))
    fout.write('{:2s} {:>12.6f} {:>12.6f} {:>12.6f}\n'.format(atomtype, x[0], x[1], x[2]))

def read_boxdata():
  with open("../warehouse/BOXDATA") as f:
    ll = [l.strip() for l in f.readlines()]
    res = [list(sub) for ele, sub in groupby(ll, key = bool) if ele] 
    k = [item[0] for item in res]
    v = [item[1:] for item in res]
    v = [[i.split(' ') if len(i.split(' '))>1 else i for i in sl] for sl in v]
    return dict(zip(k,v))

BOXDATA = read_boxdata()
trans_val = float(BOXDATA['$ZTRASL'][0].replace('d','.'))


def open_file(ff):

  f= open(ff)
  g= open(gg,'w')
  h= open(hh,'w')
  #j= open(jj,'w')
  #mass_arr, atomtypes = memoize_mass(f)   
  for i in count(1):
    try:
      xyz = read_xyz(f)
 
    except ValueError as e:
      print("DONE")
      break

    if i%1==0: print(f"x{i}")
    t = calcul_CM(xyz, kind='Si')
    water_mols = water_molecules(xyz)
    water_com = calcul_CM(water_mols, kind='water')
    #print(water_mols)
    #t_trans = translate(water_molecules(xyz),-trans_val)
    write_xyz(g, t, title="", atomtypes=xyz.atomtypes) 
    write_xyz(h, water_com, atomtypes=xyz.atomtypes) 
     
  f.close()
  g.close()
  h.close()

open_file(ff)
