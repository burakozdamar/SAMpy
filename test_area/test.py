import sys
import numpy as np

sys.path.append('/Users/burakozdamar/myGitRep/SFGpy')

from utils.xyz_utils import read_xyz

ff = "/Users/burakozdamar/brites-tmp/PEG8_traj2_pos70to130ps-10steps.xyz"

#xyz = read_xyz(ff)

with open(ff) as f:
    xyz = read_xyz(f)
    print(xyz)
