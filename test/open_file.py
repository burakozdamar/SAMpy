#!/usr/bin/env/python3

import sys
import numpy as np

sys.path.append('..')

from utils.xyz_utils import read_xyz

ff = "PEG8_traj2_pos70to130ps-10steps.xyz"
mass = {'H': 1.00794,'C1':1, 'C':12.0107, 'O':15.9994,'Si':28.0855,'Cl':35.4527,'K':39.0983,'Al':26.981539}
pbc_ = np.array([13.386,13.286,85])#*1044).reshape(-1,3)

def memoize_mass():
  print("memo")
  with open(ff) as f:
    xyz = read_xyz(f)
    atomtypes = xyz.atomtypes
    mass_arr = [mass[atom] for atom in atomtypes]
    mass_arr = np.array(mass_arr)
  return mass_arr

mass_arr = memoize_mass() 

def translate(arr1,arr2):
  return arr1-arr2

def calcul_CM(xyz):
  coords = xyz.coords
  com = np.average(coords, axis=0, weights=mass_arr)
  t = translate(coords, com)
  print("coords", coords)
  print("com", com)
  print("t",t)
  pbc(t)
  #print(com)
  #print(coords - com)

def pbc(arr):
  print(arr[0,0])
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

  print(arr[0])
  #print(arr[arr < -pbc_/2][:5])
  #print(pbc_/2)

def open_file(ff):

  with open(ff) as f:
    for i in range(100):
      try:
        xyz = read_xyz(f)
        print(f"Finished reading {i}")
      except:
        print("DONE")
        break
  
      print("calc")
      calcul_CM(xyz)

open_file(ff)
