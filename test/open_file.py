#!/usr/bin/env/python3

import sys
import numpy as np

sys.path.append('..')

ff = sys.argv[1]

from utils.xyz_utils import read_xyz
from collections import Counter 

#ff = "../warehouse/mini.xyz"
#ff = "../warehouse/PEG.xyz"
#ff = "../warehouse/test.xyz"

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
    
  print(Counter(atomtypes))

  #mass_arr = np.array([mass[atom] for atom in atomtypes if atom!="C" else 0])
  #mass_arr = np.array([mass[atom] for atom in atomtypes])
  print("orig", mass_arr)
  #mass_arr = mass_arr[np.nonzero(mass_arr)]
  print("modif", mass_arr)
  return mass_arr

mass_arr = memoize_mass(ff) 

def translate(arr1,arr2):
  return arr1-arr2

def calcul_CM(xyz):
  coords = xyz.coords
  atomtypes = xyz.atomtypes
  com = np.average(coords, axis=0, weights=mass_arr)
  t = translate(coords, com)
  print("com", com)
  print('t',t)
  #t_prime = pbc(t)
  #print("t'",t_prime)
  #return t

def pbc(arr):
  arr = np.where(arr < (-pbc_/2), arr+(pbc_/2), arr)
  arr = np.where(arr > (+pbc_/2), arr-(pbc_/2), arr)
  return arr 



def pbc1(arr):

  print("before:", arr[0])
  if arr[0,0] < -pbc_[0]/2:
      arr[0,0] += pbc_[0]
  if arr[0,0] > pbc_[0]/2:
      arr[0,0] -= pbc_[0]

  if arr[0,1] < -pbc_[1]/2:
      arr[0,1] += pbc_[1]
  if arr[0,1] > pbc_[1]/2:
      arr[0,1] -= pbc_[1]

  if arr[0,2] < -pbc_[2]/2:
      arr[0,2] += pbc_[2]
  if arr[0,2] > pbc_[2]/2:
      arr[0,2] -= pbc_[2]

  print("after",arr[0])
  #print(arr[arr < -pbc_/2][:5])
  #print(pbc_/2)

def write_to_file():
  pass

def open_file(ff):

  with open(ff) as f:
    for i in range(2):
      try:
        xyz = read_xyz(f)
        print(f"{i}")
      except:
        print("DONE")
        break
  
      calcul_CM(xyz)

open_file(ff)
