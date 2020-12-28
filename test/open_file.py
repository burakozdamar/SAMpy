#!/bin/envs/python3

import sys
import numpy as np

sys.path.append('/Users/burakozdamar/myGitRep/SFGpy')

from utils.xyz_utils import read_xyz
#import utils.xyz_utils

ff = "/Users/burakozdamar/brites-tmp/PEG8_traj2_pos70to130ps-10steps.xyz"
mass = {'H': 1.00794,'C':12.0107,'O':15.9994,'Si':28.0855,'Cl':35.4527,'K':39.0983,'Al':26.981539}

def memoize_mass():
  print("memo")
  with open(ff) as f:
    xyz = read_xyz(f)
    coords = xyz.coords
    atomtypes = xyz.atomtypes
    mass_arr = np.zeros((len(coords),))
    for idx, atom in enumerate(atomtypes):
      mass_arr[idx] = mass[atom]
  return mass_arr

mass_arr = memoize_mass() 

def calcul_CM(xyz):
  coords = xyz.coords
  atomtypes = xyz.atomtypes
  
  x_sum = np.sum(coords[:,0] * mass_arr)
  y_sum = np.sum(coords[:,1] * mass_arr)
  z_sum = np.sum(coords[:,2] * mass_arr)

  print(x_sum,y_sum,z_sum) 
  if x_sum/np.sum(mass_arr) > 13.386/2:
      print( x_sum/np.sum(mass_arr) -  13.386/2)
  print(x_sum/np.sum(mass_arr))
  print(y_sum/np.sum(mass_arr))
  print(z_sum/np.sum(mass_arr))

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
