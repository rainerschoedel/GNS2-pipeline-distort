#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul  3 16:26:43 2022

@author: alvaromartinez1

Modified by Rainer Schoedel on 3 May 2024 to use all chips simultaneously
"""

# =============================================================================
# Load all the dejitter images in the same cube for ocular inspection
# =============================================================================

import numpy as np
from scipy import ndimage
from astropy.io import fits
import glob
import sys
from astropy.io import fits
import os


common_path = '/home/data/GNS/2021/H/20HB/cubes/'
pruebas = '/home/data/GNS/2021/H/20HB/pruebas/'
basename = 'clean_GC_H_F20.H.2022-08-02T00_51_44_'
n_ims = 490 #number of exposures
i=0

n1 = 4840
n2 = 4836

nT = 0
cube = np.zeros((100,n2,n1))
for cc in range(100):
    n_exp = str(cc+1)
    imExist = os.path.exists(common_path + basename + n_exp + '.01.fits')
    if imExist:
        data_c1 = fits.getdata(common_path + basename + n_exp + '.01.fits')
        data_c2 = fits.getdata(common_path + basename + n_exp + '.02.fits')
        data_c3 = fits.getdata(common_path + basename + n_exp + '.03.fits')
        data_c4 = fits.getdata(common_path + basename + n_exp + '.04.fits')
        cube[nT,:,:] = data_c1 + data_c2 + data_c3 + data_c4
        nT = nT + 1
    print(f'Processed exposuse {n_exp}...')
fits.writeto(pruebas + 'check100.fits',cube[0:nT,:,:], overwrite=True)


nT = 0
cube = np.zeros((100,n2,n1))
for cc in range(100):
    n_exp = str(100 + cc+1)
    imExist = os.path.exists(common_path + basename + n_exp + '.01.fits')
    if imExist:
        data_c1 = fits.getdata(common_path + basename + n_exp + '.01.fits')
        data_c2 = fits.getdata(common_path + basename + n_exp + '.02.fits')
        data_c3 = fits.getdata(common_path + basename + n_exp + '.03.fits')
        data_c4 = fits.getdata(common_path + basename + n_exp + '.04.fits')
        cube[nT,:,:] = data_c1 + data_c2 + data_c3 + data_c4
        nT = nT + 1
    print(f'Processed exposuse {n_exp}...')
fits.writeto(pruebas + 'check200.fits',cube[0:nT,:,:], overwrite=True)


nT = 0
cube = np.zeros((100,n2,n1))
for cc in range(100):
    n_exp = str(200+cc+1)
    imExist = os.path.exists(common_path + basename + n_exp + '.01.fits')
    if imExist:
        data_c1 = fits.getdata(common_path + basename + n_exp + '.01.fits')
        data_c2 = fits.getdata(common_path + basename + n_exp + '.02.fits')
        data_c3 = fits.getdata(common_path + basename + n_exp + '.03.fits')
        data_c4 = fits.getdata(common_path + basename + n_exp + '.04.fits')
        cube[nT,:,:] = data_c1 + data_c2 + data_c3 + data_c4
        nT = nT + 1
    print(f'Processed exposuse {n_exp}...')
fits.writeto(pruebas + 'check300.fits',cube[0:nT,:,:], overwrite=True)


nT = 0
cube = np.zeros((100,n2,n1))
for cc in range(100):
    n_exp = str(300+cc+1)
    imExist = os.path.exists(common_path + basename + n_exp + '.01.fits')
    if imExist:
        data_c1 = fits.getdata(common_path + basename + n_exp + '.01.fits')
        data_c2 = fits.getdata(common_path + basename + n_exp + '.02.fits')
        data_c3 = fits.getdata(common_path + basename + n_exp + '.03.fits')
        data_c4 = fits.getdata(common_path + basename + n_exp + '.04.fits')
        cube[nT,:,:] = data_c1 + data_c2 + data_c3 + data_c4
        nT = nT + 1
    print(f'Processed exposuse {n_exp}...')
fits.writeto(pruebas + 'check400.fits',cube[0:nT,:,:], overwrite=True)


nT = 0
cube = np.zeros((100,n2,n1))
for cc in range(100):
    n_exp = str(400 + cc+1)
    imExist = os.path.exists(common_path + basename + n_exp + '.01.fits')
    if imExist:
        data_c1 = fits.getdata(common_path + basename + n_exp + '.01.fits')
        data_c2 = fits.getdata(common_path + basename + n_exp + '.02.fits')
        data_c3 = fits.getdata(common_path + basename + n_exp + '.03.fits')
        data_c4 = fits.getdata(common_path + basename + n_exp + '.04.fits')
        cube[nT,:,:] = data_c1 + data_c2 + data_c3 + data_c4
        nT = nT + 1
    print(f'Processed exposuse {n_exp}...')
fits.writeto(pruebas + 'check500.fits',cube[0:nT,:,:], overwrite=True)


    
    
    
    
    
    
    
