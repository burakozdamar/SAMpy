#!/usr/bin/env python3

import re
import sys
import numpy as np

sys.path.append('..')

ff = sys.argv[1]
hh = "py_pos_rebuilt.xyz"
jj = "py_pos_rebuilt_solid.xyz"

from utils.xyz_utils import read_xyz, read_boxdata
from collections import Counter, namedtuple
from itertools import cycle, count, groupby
from profiler.prof import profile_me

mass = {'H': 1.00794,'C':12.0107,'O':15.9994,'Si':28.0855,'Cl':35.4527,'K':39.0983,'Al':26.981539}
pbc_ = np.array([13.386,13.286,85])

BOXDATA = read_boxdata()
trans_val = abs(BOXDATA['$ZTRASL'])
num_oxygens = BOXDATA['$NO']
num_waters = num_oxygens*3

def memoize_mass(ff):
  print("MEMO")
  with open(ff) as f:
    xyz = read_xyz(f)
    
    if 'C' in xyz.atomtypes:
      mass_arr = np.array([mass[atom] if atom == 'Si' else 0 for atom in xyz.atomtypes])
    else:
      mass_arr = np.array([mass[atom] for atom in xyz.atomtypes])

  return mass_arr

def sep(xyz):
  data = xyz.data
  coords = xyz.coords
  atomtypes = xyz.atomtypes 

  com = np.average(coords, axis=0, weights=mass_arr[:])
  translated = translate(coords, com)
  translated = apply_pbc(translated)
  translated[:,2] -= trans_val

  waters = translated[:num_waters]
  solids = translated[num_waters:]

  waters = apply_pbc(waters)
  waters = rebuild_water(waters)
  
  solids = apply_pbc(solids) 
  final_tr = np.concatenate((waters, solids), axis=0) 

  return namedtuple("sepnt", ["coords", "data", "atomtypes"])(
        final_tr, data, atomtypes 
    )

def apply_pbc(arr):
  arr1 = np.where(arr < -pbc_/2, arr+pbc_, arr)
  arr = np.where(arr1 >  pbc_/2, arr-pbc_, arr1)
  return arr 

def translate(arr1,arr2):
  return arr1-arr2

def rebuild_water(t):

  O = t[::3]
  H1 = t[1::3]
  H2 = t[2::3]
  OH1 = H1-O  
  OH2 = H2-O  

  H1 += (OH1 < -pbc_/2)*pbc_ 
  H1 -= (OH1 >  pbc_/2)*pbc_ 
  H2 += (OH2 < -pbc_/2)*pbc_ 
  H2 -= (OH2 >  pbc_/2)*pbc_ 

  return t

#format is slow!
def write_xyz(fout, coords, title="", atomtypes=("A",)):
  fout.write("%d\n%s\n" % (coords.size / 3, title))
  for x, atomtype in zip(coords.reshape(-1, 3), cycle(atomtypes)):
    fout.write('{:2s} {:>12.6f} {:>12.6f} {:>12.6f}\n'.format(atomtype, x[0], x[1], x[2]))
    #fout.write("%s %.18g %.18g %.18g\n" % (atomtype, x[0], x[1], x[2]))

#@profile_me
def main(ff, boundary=1):

  f = open(ff)
  h = open(hh,'w')
  j = open(jj,'w')

  for i in count(1):
    try:
      xyz = read_xyz(f)
 
    except ValueError as e:
      print("DONE")
      break

    sep1 = sep(xyz)
    
    step_num = int(sep1.data[0])
    write_xyz(h, sep1.coords[:360], title=f"step wat = {step_num}", atomtypes=xyz.atomtypes) 
    write_xyz(j, sep1.coords , title=f"step rebu = {i}", atomtypes=xyz.atomtypes) 
     
  f.close()
  h.close()
  j.close()

mass_arr = memoize_mass(ff) 

main(ff)
#print(read_boxdata())
