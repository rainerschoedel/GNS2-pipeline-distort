#!/usr/bin/env python
# coding: utf-8

# In[13]:


import numpy as np
import astroalign as aa
from astropy.io import fits
from astropy.table import Table
import pandas as pd
import sys

t_exp = 3.32
ZP = 26.3

# number of brightest stars in the lists that wil be used
# (this will speed up things a lot)
n_bright = 1000


#get parameters from script call
field = int(sys.argv[1])
chip  = str(sys.argv[2])
band  = sys.argv[3]


#define paths
VVV_path  ='/home/data/VVV/Fields/J/'
#lnx_path  = f'/home/data/GNS/2021/{band}/{field}HB/ims/'
#lnx_stars = f'/home/data/GNS/2021/{band}/{field}HB/data/'
data_path = f'/home/data/GNS/2021/{band}/{field}HB/data/'
#cube_path = f'/home/data/GNS/2021/{band}/{field}HB/cubes/'

#Read VVV data
img_vvv = fits.getdata(VVV_path + f'Field{field}.fits.gz')
vvv_img = pd.DataFrame(np.array(img_vvv).byteswap().newbyteorder())
xsize = vvv_img.shape[0]
ysize = vvv_img.shape[1]
xh = round(xsize/2)
yh = round(ysize/2)

#read VVV table
#The table is sorted from bright to faint KS magnitudes
#Use only the brightest stars to save time
print(f'Size of VVV image is {xsize} x {ysize}.')
tab_vvv = Table.read(VVV_path + f'Field{field}_stars.txt',format='ascii')
select = np.nonzero((tab_vvv['Ks'] > 10) & (tab_vvv['Ks'] < 12))
tab_vvv = tab_vvv[select]
n_vvv = len(tab_vvv)
x_vvv = tab_vvv['x']
y_vvv = tab_vvv['y']
m_vvv = tab_vvv['Ks']
print(f'There are {n_vvv} stars in the VVV image.') 

#crop list for each chip
edge = 150
if (chip == '1'):
    idx = np.nonzero((x_vvv < xh) & (y_vvv < yh) & (x_vvv > edge) & (y_vvv > edge))
elif (chip == '2'):
    idx = np.nonzero((x_vvv > xh) & (y_vvv < yh) & (x_vvv < (xsize-edge))& (y_vvv > edge))
elif (chip == '3'):
    idx = np.nonzero((x_vvv < xh) & (y_vvv > yh) & (x_vvv > edge) &  (y_vvv < (ysize-edge)))
elif (chip == '4'):
    idx = np.nonzero((x_vvv > xh) & (y_vvv > yh) & (x_vvv < (xsize-edge)) & (y_vvv < (ysize-edge)))
x_vvv = x_vvv[idx]
y_vvv = y_vvv[idx]
m_vvv = m_vvv[idx]
tab_vvv = Table([x_vvv,y_vvv,m_vvv],names=['x','y','m'])
tab_vvv.write('tmp/vvv_prelim.txt', format='ascii',overwrite=True)
#The VVV list can contain the same star with slightly different values
#eliminate these repetitions
tol = 1
sel = np.ones(len(x_vvv))
for s in range(len(x_vvv)):
    for h in range((s+1),len(x_vvv)):
        d = np.sqrt(((x_vvv[s])-x_vvv[h])**2 + ((y_vvv[s])-y_vvv[h])**2)
        if (d < tol):
            sel[s] = 0
            print(f'Deselected star {s+1} and {h+1}.')
sel = sel.astype('bool')
x_vvv = x_vvv[sel]
y_vvv = y_vvv[sel]
m_vvv = m_vvv[sel]

print(f'There are {len(x_vvv)} stars in the corresponding section of the VVV image.')
print(f'x-range of stars in VVV image from {np.min(x_vvv)} to {np.max(x_vvv)}.')
print(f'y-range of stars in VVV image from {np.min(y_vvv)} to {np.max(y_vvv)}.')
vvv_xy = np.array([[x, y] for x, y in zip(x_vvv, y_vvv)], dtype="float64")
#print(vvv_xy)
tab_vvv = Table([x_vvv,y_vvv,m_vvv],names=['x','y','m'])
tab_vvv.write('tmp/vvv.txt', format='ascii',overwrite=True)
regfile = open(VVV_path + 'align.reg','w')
regfile.write('# Region file format: DS9 version 4.1 \n')
regfile.write('global color=white dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1 \n')
regfile.write('physical \n')
for i in range(len(x_vvv)):
    line = f'circle({round(x_vvv[i])},{round(y_vvv[i])},10) \n'
    regfile.write(line)
regfile.close()

#Read HAWK-I data
#limit number of stars to the brightest n_bright
#The input list is supposed to be ordered by descending brightness
tab_img = Table.read(data_path + 'stars_' + chip + '_holo.txt',format='ascii')
m_img = ZP - 2.5*np.log10(tab_img['f']/t_exp)
select = np.nonzero((m_img < 14))
x_img = tab_img['x'][select]
y_img = tab_img['y'][select]
#crop list at the edges
edge = 150
xhi = np.max(x_img) - edge
yhi = np.max(y_img) - edge
idx = np.nonzero((x_img < xhi) & (y_img < yhi) & (x_img > edge) & (y_img > edge))
x_img = x_img[idx]
y_img = y_img[idx]
img_xy = np.array([[x, y] for x, y in zip(x_img, y_img)], dtype="float64")
print(f'There are {len(x_img)} stars in the corresponding HAWK-I image.')
print(f'x-range of stars in HAWK-I image from {np.min(x_img)} to {np.max(x_img)}.')
print(f'y-range of stars in HAWK-I image from {np.min(y_img)} to {np.max(y_img)}.')
tab_bright = Table([x_img,y_img],names=['x','y'])
tab_bright.write('tmp/hawki.txt', format='ascii',overwrite=True)
regfile = open(data_path + 'align.reg','w')
regfile.write('# Region file format: DS9 version 4.1 \n')
regfile.write('global color=white dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1 \n')
regfile.write('physical \n')
for i in range(len(x_img)):
    line = f'circle({round(x_img[i])},{round(y_img[i])},10) \n'
    regfile.write(line)
regfile.close()

#TRANSFORMATION INTO VVV
#save ordered lists of common stars
#===================================
#Determine transformation into VVV orig
p, (pos_img, pos_img_t) = aa.find_transform(img_xy, vvv_xy, max_control_points=200)
print("\nTranformation matrix:\n{}".format(p.params))
np.savetxt(data_path + f'aa_stars_hawki_'+chip+'.txt',pos_img)
np.savetxt(data_path + f'aa_stars_vvv_'+chip+'.txt',pos_img_t)
