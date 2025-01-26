
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 20 June 2024

@author: Rainer Schoedel


Notes: I do not use the glob.glob function to obtain the paths to the files, because the 
latter will return files in random order and because it is also more complicated to match images and weights.

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
chip = int(sys.argv[2])

in_path = '/home/data/GNS/2021/H/20HB/cubes/'
wpath = '/home/data/GNS/2021/H/20HB/weights/'
out_path = '/home/data/GNS/2021/H/20HB/subcubes/'
basename = 'clean_GC_H_F20.H.2022-08-02T00_51_44_'

#read global minimum sky value to corrrect images for negativities
with open(in_path+'sky_offset.txt') as file:
    sky_min = float(file.read())

n_ims_max = 500 #maximum number of exposures

x_cube = 600 # xaxis length of sub-cube
y_cube = 600 # yaxis length of sub-cube

#The x_large and y_large values will usually NOT be the same
#Choose a suitably large value for both of them 
#This may require padding images with zeros and/or may also imply minor loss near the edges.
#Here I choose 4800 (the actual sizes are x: 4840 and y 4836)

x_large = 4800  # xaxis length of large cube
y_large = 4800  # yaxis length of large cube
valid_frac = 0.4 # minimum fraction of valid pixels required

x_sub_shift = int(x_cube/2)
y_sub_shift = int(y_cube/2)
nx = int(x_large//x_sub_shift - 1)
ny = int(y_large//y_sub_shift - 1)

n_grid = nx * ny
npix_sub = x_cube * y_cube

#First clean target directories to avoid appending data to already existing files
for i_x in range(nx):
    for i_y in range(ny):

        outnam = out_path + 'chip' + str(chip) + '/cube_' + str(i_x) + '_' + str(i_y) + '.fits'
        outmasknam = out_path + 'chip' + str(chip) + '/masks_' + str(i_x) + '_' + str(i_y) + '.fits'
        if os.path.exists(outnam):
            os.remove(outnam)
        if os.path.exists(outmasknam):
            os.remove(outmasknam)

#Loop over exposures
#count = 0
for i_f in range(n_ims_max):

    n_exp = str(i_f + 1)
    impath = in_path + basename + n_exp + '.0' + str(chip) + '.fits'
    maskpath = wpath + basename + n_exp + '.0' + str(chip) + '.weight.fits'
    imExist = os.path.exists(impath)
    mExist = os.path.exists(maskpath)

    if (imExist and mExist):

        print(f'Now working on image {impath}.')
        im = fits.getdata(impath)
        hdr = fits.getheader(impath)
        mm = fits.getdata(maskpath)
        mm = (mm > 4.e17)
        mask = mm.astype(int)
        im = im - sky_min
        im = im * mask
        im = im * hdr['GAIN']

        #Loop over sub-fields
        for i_x in range(nx):
            for i_y in range(ny):

                outnam = out_path + 'chip' + str(chip) + '/cube_' + str(i_x) + '_' + str(i_y) + '.fits'
                outmasknam = out_path + 'chip' + str(chip) + '/masks_' + str(i_x) + '_' + str(i_y) + '.fits'

                xlo = i_x * x_sub_shift
                xhi = xlo + x_cube
                ylo = i_y * y_sub_shift
                yhi = ylo + y_cube

                # Field valid or completely/largely masked?
                out_im = im[ylo:yhi,xlo:xhi]
                out_mask = mask[ylo:yhi,xlo:xhi]
                wt = np.sum(out_mask)
                if ((wt/npix_sub) > valid_frac):
#                    print(f'Subcube {i_x}, {i_y}:')
#                    print(f'Mask sum: {wt}, sub sum {npix_sub}, and fraction {wt/npix_sub}.')
#                    print('*****')
#                    fits.writeto('tmp/sub'+str(count)+'.fits',out_im, overwrite=True)
#                    fits.writeto('tmp/submask'+str(count)+'.fits',out_mask, overwrite=True)
                    fits.append(outnam, out_im)
                    fits.append(outmasknam, out_mask)

#Convert MEFs to FITS cubes
for i_x in range(nx):
    for i_y in range(ny):

        outnam = out_path + 'chip' + str(chip) + '/cube_' + str(i_x) + '_' + str(i_y) + '.fits'
        outmasknam = out_path + 'chip' + str(chip) + '/masks_' + str(i_x) + '_' + str(i_y) + '.fits'

        if os.path.exists(outnam):
            with fits.open(outnam) as hdul:
                n_ext = len(hdul)
                cube = np.zeros((n_ext,y_cube,x_cube))
                l = 0
                for hdr in hdul:
                    cube[l,:,:] = hdul[l].data
                    l += 1
                fits.writeto(outnam,cube,overwrite=True)

        if os.path.exists(outmasknam):
            with fits.open(outmasknam) as hdul:
                n_ext = len(hdul)
                cube = np.zeros((n_ext,y_cube,x_cube))
                l = 0
                for hdr in hdul:
                    cube[l,:,:] = hdul[l].data
                    l += 1
            fits.writeto(outmasknam,cube,overwrite=True)

print(f'Finished chip {chip}.')    
    
print(f'Sky minimum is {sky_min}.')
