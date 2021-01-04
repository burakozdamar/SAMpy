#!/usr/bin/env/python3

import re
import sys
import numpy as np

sys.path.append('..')

ff = sys.argv[1]
gg = "py_pos_recentered.xyz"
hh = "py_pos_rebuilt.xyz"
jj = "py_pos_rebuilt_solid.xyz"

from utils.xyz_utils import read_xyz, read_boxdata
from collections import Counter, namedtuple
from itertools import cycle, count, groupby

mass = {'H': 1.00794,'x':0,'C':12.0107,'O':15.9994,'Si':28.0855,'Cl':35.4527,'K':39.0983,'Al':26.981539}
pbc_ = np.array([13.386,13.286,85])

BOXDATA = read_boxdata()
trans_val = abs(float(BOXDATA['$ZTRASL'][0].replace('d','.')))

def memoize_mass(xyz):
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

def sep(xyz):
  data = xyz.data
  coords = xyz.coords
  atomtypes = xyz.atomtypes 

  if 'C' in xyz.atomtypes:
    mass_arr = np.array([mass[atom] if atom == "Si" else 0 for atom in xyz.atomtypes])  
  else:
    mass_arr = np.array([mass[atom] for atom in xyz.atomtypes])

  com = np.average(coords, axis=0, weights=mass_arr[:])
  #print(com)
  translated = translate(coords, com)
  translated = pbc(translated)

  waters = translated[:360]
  solids = translated[360:]

  waters[:,2] -= trans_val
  waters = pbc(waters)
  waters = manual_wat_pbc(waters)
  solids[:,2] -= trans_val
  solids = pbc(solids) 
  final_tr = np.concatenate((waters, solids), axis=0) 

  return namedtuple("sepnt", ["coords", "data", "atomtypes"])(
        final_tr, data, atomtypes 
    )


def pbc(arr):
 # arr = np.where(arr < -pbc_/2, arr+pbc_, arr)
 # arr = np.where(arr >  pbc_/2, arr-pbc_, arr)
 # return arr 
  arr1 = np.where(arr < -pbc_/2, arr+pbc_, arr)
  arr = np.where(arr1 >  pbc_/2, arr-pbc_, arr1)
  return arr 

def translate(arr1,arr2):
  return arr1-arr2

def translate_to_CM(xyz, kind='Si', linear_trans=0):
  data = xyz.data
  coords = xyz.coords
  atomtypes = xyz.atomtypes
  n_oxygens = BOXDATA['$NO']
  solids = coords[int(n_oxygens[0])*3:]

  if 'C' in xyz.atomtypes:
    mass_arr = np.array([mass[atom] if atom == kind else 0 for atom in xyz.atomtypes])  
  else:
    mass_arr = np.array([mass[atom] for atom in xyz.atomtypes])

  com = np.average(coords, axis=0, weights=mass_arr[:coords.shape[0]])
  translated = translate(coords, com)

  #water_mols = water_molecules(translated, translate=trans_val, rebuilt=True)
  solids -= linear_trans
  
  t_solids = pbc(solids)

  #final_tr = np.concatenate((water_mols.coords, t_solids), axis=0) 
  final_tr = 0
  #translated = pbc(translated)

  return namedtuple("trans_to_CM", ["coords", "data", "atomtypes"])(
        final_tr, data, atomtypes 
    )
  #return translated

def water_molecules(xyz, translate=0, rebuilt=False):
  data = xyz.data
  atomtypes = xyz.atomtypes
  coords = xyz.coords
  n_oxygens = BOXDATA['$NO']
  coords = coords[:int(n_oxygens[0])*3]  
  atomtypes = atomtypes[:int(n_oxygens[0])*3]  
  coords[:,2] -= translate

  if rebuilt:
    #print("rebuilt", coords[0])
    coords = pbc(coords)
    #print(coords[0])
    #coords = water_pbc(coords)
    #coords = test_pbc(coords)
    #print(coords[0])
  #print(coords.shape)

  return namedtuple("waters", ["coords", "data", "atomtypes"])(
        coords, data, atomtypes 
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


def water_pbc(t):
  print(t[:10])
  O = t[::3]
  H1 = t[1::3]
  H2 = t[2::3]
  OH1 = H1-O  
  OH2 = H2-O  
  print(OH1[:10])

  arr = np.where(OH1 < -pbc_/2, H1+pbc_, H1)
  print(t[:3])
  arr = np.where(OH1 >  pbc_/2, arr-pbc_, H1)
  
  print(arr[:10])
  arr = np.where(OH2 < -pbc_/2, arr+pbc_, H2)
  arr = np.where(OH2 >  pbc_/2, arr-pbc_, H2)
  #not sure if it works
  #print(t[10])
  return t

def manual_wat_pbc(t):
  print('t',t[:10])
  O = t[::3]
  H1 = t[1::3]
  H2 = t[2::3]
  OH1 = H1-O  
  OH2 = H2-O  
  
  for i in range(OH1.shape[0]):
    for j in range(OH1.shape[1]):
      if OH1[i][j] < -pbc_[j]/2:
        H1[i][j]+pbc_[j]
      elif OH1[i][j] > +pbc_[j]/2:
        H1[i][j]-pbc_[j]

      if OH2[i][j] < -pbc_[j]/2:
        H2[i][j]+pbc_[j]
      elif OH2[i][j] > +pbc_[j]/2:
        H2[i][j]-pbc_[j]

  #H1[0][0]=123456789
  return t

def test_pbc(t):
  O = t[::3]
  H1 = t[1::3]
  H2 = t[2::3]
  
  O[0,0]=9999

  OH1 = pbc(O-H1)  
  OH2 = pbc(O-H2)  
  print(t)
  return t  
  
def write_xyz(fout, coords, title="", atomtypes=("A",)):
  fout.write("%d\n%s\n" % (coords.size / 3, title))
  for x, atomtype in zip(coords.reshape(-1, 3), cycle(atomtypes)):
    #fout.write("%s %.8g %.8g %.8g\n" % (atomtype, x[0], x[1], x[2]))
    fout.write('{:2s} {:>12.6f} {:>12.6f} {:>12.6f}\n'.format(atomtype, x[0], x[1], x[2]))


def main(ff, boundary=1):

  f = open(ff)
  g = open(gg,'w')
  h = open(hh,'w')
  j= open(jj,'w')
  #mass_arr, atomtypes = memoize_mass(f)   
  for i in count(1):
    try:
      xyz = read_xyz(f)
 
    except ValueError as e:
      print("DONE")
      break

    #if i%1==0: print(f"x{i}")
    print(i)
    sep1 = sep(xyz)
    #t = translate_to_CM(xyz, kind='Si', linear_trans = trans_val)
    #water_mols = water_molecules(t, translate=trans_val, rebuilt=True)
    
    step_num = int(sep1.data[0])
    #write_xyz(g, t.coords, title='x', atomtypes=xyz.atomtypes) 
    write_xyz(h, sep1.coords[:360], title=f"step = {step_num}", atomtypes=xyz.atomtypes) 
    write_xyz(j, sep1.coords , title=f"step rebu = {i}", atomtypes=xyz.atomtypes) 
     
  f.close()
  g.close()
  h.close()
  j.close()
main(ff)
