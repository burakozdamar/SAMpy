#!/usr/bin/env/python3

import sys
import numpy as np

sys.path.append('..')

from utils.xyz_utils import read_xyz
from collections import Counter 

ff = "../warehouse/PEG.xyz"
#ff = "../warehouse/test.xyz"

mass = {'H': 1.00794,'C':0,'x':12.0107,'O':15.9994,'Si':28.0855,'Cl':35.4527,'K':39.0983,'Al':26.981539}
pbc_ = np.array([13.386,13.286,85])

def memoize_mass(ff):
  print("memo")
  with open(ff) as f:
    xyz = read_xyz(f)
    atomtypes = xyz.atomtypes
    #for a in atomtypes:
    #  if a == "C":
        
    mass_arr = [mass[atom] for atom in atomtypes]
    mass_arr = np.array(mass_arr)
  return mass_arr

mass_arr = memoize_mass(ff) 

def translate(arr1,arr2):
  return arr1-arr2

def calcul_CM(xyz):
  coords = xyz.coords
  atomtypes = xyz.atomtypes
  print(Counter(atomtypes))
  print(atomtypes)
  print(mass_arr, sum(mass_arr))
  com = np.true_divide(mass_arr.sum(0),(mass_arr!=0).sum(0))
  #com = np.average(coords, axis=0, weights=mass_arr)
  t = translate(coords, com)
  #print("coords", coords)
  print("com", com)
  #print("t",t)
  pbc(t)
  #print(com)
  #print(coords - com)

def pbc(arr):
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

def open_file(ff):

  with open(ff) as f:
    for i in range(2):
      try:
        xyz = read_xyz(f)
        print(f"Finished reading {i}")
      except:
        print("DONE")
        break
  
      print("calc")
      calcul_CM(xyz)

open_file(ff)
