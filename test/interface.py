#!/usr/bin/env/python3

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

BOXDATA = read_boxdata()
pbc_ = np.array([13.386,13.286,85])
trans_val = abs(float(BOXDATA['$ZTRASL'][0].replace('d','.')))
num_oxygens = int(BOXDATA['$NO'][0])
num_waters = int(BOXDATA['$NO'][0])*3

def memoize_mass(ff):
  print("MEMO")
  with open(ff) as f:
    xyz = read_xyz(f)
    
    if 'C' in xyz.atomtypes:
      mass_arr = np.array([mass[atom] if atom == 'Si' else 0 for atom in xyz.atomtypes])
    else:
      mass_arr = np.array([mass[atom] for atom in xyz.atomtypes])

  return mass_arr

def process(xyz):
  data = xyz.data
  coords = xyz.coords
  atomtypes = xyz.atomtypes 
  n_wat = int(BOXDATA['$NO'][0])*3
  x,y,z = grid(pbc_)
  print(x)
   

 # return namedtuple("sepnt", ["coords", "data", "atomtypes"])(
 #       final_tr, data, atomtypes 
 #   )

def apply_pbc(arr):
  arr1 = np.where(arr < -pbc_/2, arr+pbc_, arr)
  arr = np.where(arr1 >  pbc_/2, arr-pbc_, arr1)
  return arr 

def translate(arr1,arr2):
  return arr1-arr2

#format is slow!
def write_xyz(fout, coords, title="", atomtypes=("A",)):
  fout.write("%d\n%s\n" % (coords.size / 3, title))
  for x, atomtype in zip(coords.reshape(-1, 3), cycle(atomtypes)):
    fout.write('{:2s} {:>12.6f} {:>12.6f} {:>12.6f}\n'.format(atomtype, x[0], x[1], x[2]))
    #fout.write("%s %.18g %.18g %.18g\n" % (atomtype, x[0], x[1], x[2]))

def grid(pbc):
  E = 2.4
  d = np.array([0.5, 0.5, 0.25])
  n = np.array(list(map(round, (pbc/d)[:])))
  f = - (pbc - d)/2
  k = (n-1)*d+f

  x = np.arange(f[0], k[0]+0.5 , 0.5)
  y = np.arange(f[1], k[1]+0.5 , 0.5)
  z = np.arange(f[2], k[2]+0.25 , 0.25)

  #print (d, n, f, (n-1)*d+f)
  return x,y,z

 # r = np.sqrt(xdiff**2+ydiff**2+zdiff**2)
 # if (r<=3*E):
 #   p = p + exp(-r**2/(2*E**2))/((2*np.pi*E**2)**1.5)




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

    proc = process(xyz)
    #proc     
    #interface.xyz to write
    #grid_interface to write
    #write_xyz(h, sep1.coords[:360], title=f"step wat = {step_num}", atomtypes=xyz.atomtypes) 
    #write_xyz(j, sep1.coords , title=f"step rebu = {i}", atomtypes=xyz.atomtypes) 
     
  f.close()
  h.close()
  j.close()

mass_arr = memoize_mass(ff) 

main(ff)
#print(read_boxdata())
