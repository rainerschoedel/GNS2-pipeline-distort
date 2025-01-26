
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
import glob
import os
import sys

field = int(sys.argv[1])

in_path = '/home/data/GNS/2021/H/20HB/cubes/'
wpath = '/home/data/GNS/2021/H/20HB/weights/'
out_path = '/home/data/GNS/2021/H/20HB/ims/'
basename = 'clean_GC_H_F20.H.2022-08-02T00_51_44_'
n_ims_max = 500 #maximum number of exposures
x_large = 4800  # xaxis length of large cube
y_large = 4800  # yaxis length of large cube

#read global minimum sky value to corrrect images for negativities
with open(in_path+'sky_offset.txt') as file:
    sky_min = float(file.read())

n_ims_max = 500

#Make mean image
#================

lxp = np.zeros((x_large,y_large))
wt = np.zeros((x_large,y_large))


#Use header from first exposure for all images
common_hdr = fits.getheader(in_path + basename + '1.01.fits')

#Loop over exposures
for i_f in range(n_ims_max):
    for i_c in range(4):
        n_exp = str(i_f + 1)
        impath = in_path + basename + n_exp + '.0' + str(i_c+1) + '.fits'
        maskpath = wpath + basename + n_exp + '.0' + str(i_c+1) + '.weight.fits'
        imExist = os.path.exists(impath)
        mExist = os.path.exists(maskpath)

        if (imExist and mExist):

            print(f'Now working on image {impath}.')
            im = fits.getdata(impath)
            hdr = fits.getheader(impath)
            im = im - sky_min
            im = im * hdr['GAIN']
            mm = fits.getdata(maskpath)
            mm = (mm > 4.e17)
            mask = mm.astype(int)
            im = im[:y_large,:x_large]
            mask = mask[:y_large,:x_large]
            lxp = lxp + im*mask
            wt = wt + mask

good = np.nonzero(wt > 0)
lxp[good] = lxp[good]/wt[good]
fits.writeto(out_path + 'lxp.fits.gz', lxp, header=common_hdr, overwrite=True)
fits.writeto(out_path + 'lxp_wt.fits.gz', wt, header=common_hdr, overwrite=True)
print(f'Finished long exposure for chip {i_c+1}.')    

#Make std image
#================

lxp = fits.getdata(out_path + 'lxp.fits.gz')
noise = np.zeros((x_large,y_large))
wt = np.zeros((x_large,y_large))

#Loop over exposures
for i_f in range(n_ims_max):
    for i_c in range(4):

        n_exp = str(i_f + 1)
        impath = in_path + basename + n_exp + '.0' + str(i_c+1) + '.fits'
        maskpath = wpath + basename + n_exp + '.0' + str(i_c+1) + '.weight.fits'
        imExist = os.path.exists(impath)
        mExist = os.path.exists(maskpath)

        if (imExist and mExist):

            print(f'Now working on image {impath}.')
            im = fits.getdata(impath)
            hdr = fits.getheader(impath)
            im = im - sky_min
            im = im * hdr['GAIN']
            mm = fits.getdata(maskpath)
            mm = (mm > 4.e17)
            mask = mm.astype(int)
            im = im[:y_large,:x_large]
            mask = mask[:y_large,:x_large]
            good = np.nonzero(mask > 0)
            noise[good] = noise[good] + (im[good] - lxp[good])**2
            wt = wt + mask

good = np.nonzero(wt > 0)
noise[good] = np.sqrt(noise[good]/(wt[good]**2))
fits.writeto(out_path + 'lxp_sigma.fits.gz', noise, header=common_header, overwrite=True)
print(f'Finished noise of long exposure.')   

    
    
    
