
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 16 Jan 2025

@author: Rainer Schoedel

"""

# =============================================================================
# Read in the images provided by Hervé (after deselection of bad frames)
# and extract subcubes 
# =============================================================================

import numpy as np
from astropy.io import fits
import shutil
import glob
import sys
import os

field = int(sys.argv[1])

in_path = '/home/data/GNS/2021/H/20HB/cubes/'
filename = in_path + 'list.txt'
with open(filename) as file:
    names = [line.rstrip() for line in file]
n_ims = len(names)

minvals = list()
#Loop over exposures
#for i in range(10):
for i in range(n_ims):

    impath = in_path + names[i]
    im = fits.getdata(impath)
    v_min = np.min(im)
    minvals.append(v_min)
    print(f'Minimum value in {names[i]} is {v_min}.')

vals = np.array(minvals)
glob_min = np.min(vals)
print(f'The global minimum value is {glob_min}.')
glob_med = np.median(vals)
print(f'The median minimum value is {glob_med}.')

outfile = open(in_path + 'sky_offset.txt', 'w')
outfile.write(str(glob_med))
outfile.close()
