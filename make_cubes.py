#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul  3 16:26:43 2022

@author: Rainer Schoedel

use input from subcubes,py
convert fits files with extensions to actual fits cubes
"""
1
# =============================================================================
# Load all the dejitter images in the same cube for ocular inspection
# =============================================================================

import numpy as np
from astropy.io import fits
import glob
import sys
import os


in_path = '/home/data/GNS/2021/H/20HB/pruebas/'
out_path = '/home/data/GNS/2021/H/20HB/subcubes/'

nx = 14
ny = 14

#x and y axis size of subcubes
n1 = 600
n2 = 600

#Loop over chips
for i_c in range(4):

    #Loop over sub-fields
    for i_x in range(nx):
        for i_y in range(ny):

            #Names of files. Do they exist?
            name = 'cube_' + str(i_x) + '_' + str(i_y) + '.fits'
            imExist = os.path.exists(in_path + 'chip' + str(i_c+1) + '/' + name)
            maskname = 'holomasks_' + str(i_x) + '_' + str(i_y) + '.fits'
            maskExist = os.path.exists(in_path + 'chip' + str(i_c+1) + '/' + maskname)

            if (imExist and imExist):

                #work on data cubes
                print(f'Working on chip{i_c+1}, ' + name + '...')
                hdu = fits.open(in_path + 'chip' + str(i_c+1) + '/' + name)
                n3 = len(hdu)
                print(f'Size of cube: {n3}')
                cube = np.zeros((n3,n2,n1))
                for i in range(n3):
                    cube[i,:,:] = hdu[i].data
                fits.writeto(out_path + 'chip' + str(i_c+1) + '/'  + name, cube)

                #Now work on mask cubes
                name = 'holomasks_' + str(i_x) + '_' + str(i_y) + '.fits'
                hdu = fits.open(in_path + 'chip' + str(i_c+1) + '/' + maskname)
                n3 = len(hdu)
                print(f'Size of cube: {n3}')
                cube = np.zeros((n3,n2,n1))
                for i in range(n3):
                    cube[i,:,:] = hdu[i].data
                fits.writeto(out_path+ 'chip' + str(i_c+1) + '/'  + maskname, cube)

                with open (out_path + 'chip' + str(i_c+1) + '/list_' + str(i_x) + '_' + str(i_y) + '.txt', 'w') as f:
                    f.write('/cube_' + str(i_x) + '_' + str(i_y) + '.fits')
                with open (out_path + 'chip' + str(i_c+1) + '/masklist_' + str(i_x) + '_' + str(i_y) + '.txt', 'w') as f:
                    f.write('/masks_' + str(i_x) + '_' + str(i_y) + '.fits')

    
    
    
