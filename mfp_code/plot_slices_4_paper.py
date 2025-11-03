import numpy as np
import matplotlib.pyplot as plt
import math
from scipy import interpolate
from scipy import integrate
#import seaborn as sns
from mpl_toolkits.axes_grid1 import make_axes_locatable
from astropy import constants as const
import random
import scipy
import argparse
import pandas as pd
import random
import matplotlib.colors as colors

#global constants
MPC        = 3.08568025e+22  # m
H0         = 1.0e5/MPC       # s-1
kPC        = const.kpc.value #m
GRAVITY    = const.G.value   #m3 / (kg s2)
PROTONMASS = const.m_p.value #kg

def number_density_H(h100,overdensity, omegab, XHfrac, redshift):
    #Get nH(z) at the critical density in proper cgs units
    
    rhoc       = 3.0 * (H0*h100)**2/(8.0 * np.pi * GRAVITY) #kg m-3
    
    n_H     = rhoc*omegab*(overdensity)*XHfrac*(1.0 + redshift)**3.0 / PROTONMASS #m-3
    
    return n_H

z = ['10.0', '8.0', '7.0', '6.0', '5.4']
# create an Empty DataFrame object
df = pd.DataFrame()
for ii in z:
	file = '/home/ppxjf3/planck1_40_2048_RTzrfit_los/los2048_n5000_z'+ ii + '00.dat'
	"""_
	------------------------------------------------------------------------------------------------------------------------
	Read in binary file data 
	------------------------------------------------------------------------------------------------------------------------
	"""
	# Open the binary file
	readdata = open(file,"rb")
		
	# Header data
	ztime  = np.fromfile(readdata,dtype=np.double,count=1) # redshift
	omegaM = np.fromfile(readdata,dtype=np.double,count=1) # Omega_m (matter density)
	omegaL = np.fromfile(readdata,dtype=np.double,count=1) # Omega_L (Lambda density)
	omegab = np.fromfile(readdata,dtype=np.double,count=1) # Omega_b (baryon density)
	h100   = np.fromfile(readdata,dtype=np.double,count=1) # Hubble constant, H0 / 100 km/s/Mpc
	box100 = np.fromfile(readdata,dtype=np.double,count=1) # Box size in comoving kpc/h
	Xh     = np.fromfile(readdata,dtype=np.double,count=1) # Hydrogen fraction by mass
	nbins  = np.fromfile(readdata,dtype=np.int32,count=1)  # Number of pixels in each line of sight
	numlos = np.fromfile(readdata,dtype=np.int32,count=1)  # Number of lines of sight
	print(numlos)
	print(nbins)
	# Line of sight locations in box 
	iaxis  = np.fromfile(readdata,dtype=np.int32,count=numlos[0])  # projection axis, x=1, y=2, z=3
	xaxis  = np.fromfile(readdata,dtype=np.double,count=numlos[0]) # x-coordinate in comoving kpc/h
	yaxis  = np.fromfile(readdata,dtype=np.double,count=numlos[0]) # y-coordinate in comoving kpc/h
	zaxis  = np.fromfile(readdata,dtype=np.double,count=numlos[0]) # z-coordinate in comoving kpc/h

	# Line of sight scale
	posaxis = np.fromfile(readdata,dtype=np.double,count=nbins[0]) # comoving kpc/h
	velaxis = np.fromfile(readdata,dtype=np.double,count=nbins[0]) # km/s

	# Gas density, rho/<rho>
	density = np.fromfile(readdata,dtype=np.double,count=nbins[0]*numlos[0])

	# H1 fraction, fH1 = nH1/nH
	H1frac  = np.fromfile(readdata,dtype=np.double,count=nbins[0]*numlos[0])

	# Temperature, K
	temp_1   = np.fromfile(readdata,dtype=np.double,count=nbins[0]*numlos[0])

	# Peculiar velocity, km/s
	vpec    = np.fromfile(readdata,dtype=np.double,count=nbins[0]*numlos[0])

	# Close the binary file
	readdata.close()
 
	n = number_density_H(h100,density, omegab, Xh, ztime)*1e-6 #cm-3
	n_HI = n*H1frac #cm-3
	if ii == '10.0':
		df['10.0'] = np.log10(n_HI)
	if ii == '8.0':
		df['8.0'] = n_HI
	if ii == '7.0':
		df['7.0'] = n_HI
	if ii == '6.0':
		df['6.0'] = n_HI
	if ii == '5.4':
		df['5.4'] = n_HI
  
"""  
im3 = ax[1,1].imshow(NHI_slice, extent=[0,float(box100),float(box100),0], cmap='viridis')
cbar3 = fig.colorbar(im3,ax=ax[1,1])
ax[1,1].set_ylabel('z [cMpc/h]')
ax[1,1].set_xlabel('y [cMpc/h]')
cbar3.set_label(r'log$_{10}[n_{H1}/n_H ]$')
ax[1,1].invert_yaxis()"""