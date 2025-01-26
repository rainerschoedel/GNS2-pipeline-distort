import numpy as np
from astropy.io import fits

in_path = '/home/data/GNS/2021/H/20HB/ims/'
out_path = '/home/data/GNS/2021/H/20HB/data/'

n1 = 2700
n2 = 2700
n1_large = 4800
n2_large = 4800

# chip 1
im = fits.getdata(in_path + 'lxp_chip1.fits.gz')    
noise = fits.getdata(in_path + 'lxp_chip1_sigma.fits.gz')
im = im[:n2,:n1]
noise = noise[:n2,:n1]
fits.writeto(out_path + 'lxp_chip1.fits.gz', im, overwrite=True)
fits.writeto(out_path + 'lxp_chip1_sigma.fits.gz', noise, overwrite=True)
print('Finished chip1.')

# chip 2
im = fits.getdata(in_path + 'lxp_chip2.fits.gz')    
noise = fits.getdata(in_path + 'lxp_chip2_sigma.fits.gz')
im = im[:n2,(n1_large-n1):]
noise = noise[:n2,(n1_large-n1):]
fits.writeto(out_path + 'lxp_chip2.fits.gz', im, overwrite=True)
fits.writeto(out_path + 'lxp_chip2_sigma.fits.gz', noise, overwrite=True)
print('Finished chip 2.')

# chip 3
im = fits.getdata(in_path + 'lxp_chip3.fits.gz')    
noise = fits.getdata(in_path + 'lxp_chip3_sigma.fits.gz')
im = im[(n2_large-n2):,:n1]
noise = noise[(n2_large-n2):,:n1]
fits.writeto(out_path + 'lxp_chip3.fits.gz', im, overwrite=True)
fits.writeto(out_path + 'lxp_chip3_sigma.fits.gz', noise, overwrite=True)
print('Finished chip 3.')

# chip 4
im = fits.getdata(in_path + 'lxp_chip4.fits.gz')    
noise = fits.getdata(in_path + 'lxp_chip4_sigma.fits.gz')
im = im[(n2_large-n2):,(n1_large-n1):]
noise = noise[(n2_large-n2):,(n1_large-n1):]
fits.writeto(out_path + 'lxp_chip4.fits.gz', im, overwrite=True)
fits.writeto(out_path + 'lxp_chip4_sigma.fits.gz', noise, overwrite=True)
print('Finished chip 4.')
