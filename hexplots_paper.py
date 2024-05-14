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
import seaborn as sns

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

def find_nearest(array, value):
		array = np.asarray(array)
		idx = (np.abs(array - value)).argmin()
		return idx  

def pdf_calc (bin_min, bin_max, nbins, data, dr):
	step = (bin_max - bin_min)/nbins
	#print(step)
	bins = np.arange(bin_min, bin_max, step)
	N = np.array([])
	for j in range (0, nbins-1):
		
		Ni =np.where((data>= bins[j]) & (data < bins[j+1]))
		N = np.append(N,len(Ni[0]))
	
	N_tot = np.sum(N)
	#print(step)
	pdf = N/(step*dr)
	centers=0.5*(bins[:-1]+bins[1:])
	return pdf, bins, centers

def normalisation_factor(v, itter, h100):
    R = v/(H0*h100) #km
    R_m = R*100 #m
    R_MPc = R_m/MPC
    return R_MPc*itter

def trap_integ(x, y):
    if len(x) == 0:
        return 0
    idx = np.argsort(x)
    x = x[idx]
    y = y[idx]
    sum = 0 
    for i in range(1,len(y)-2):
        sum += ((y[i-1] + y[i])/2 )*(x[i]- x[i-1])
    return sum

def logmean(x):
    return np.mean(np.log10(x))

def loglen(x):
    return np.log10(len(x))

window = 50.0    
z = ['10.0', '6.0', '5.4']
N_HI_100 = np.array([])
N_HI_60 = np.array([])
N_HI_54 = np.array([])
for ii in z:
	file = '/home/ppxjf3/planck1_40_2048_RTzrfit_los/los2048_n5000_z'+ ii + '00.dat'
	"""_
	------------------------------------------------------------------------------------------------------------------------
	Read ilues over which the contour is drawn. Color-mapping is controlled by cmap, norm, vmin, and vmax.n binary file data 
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

	""""
	---------------------------------------------------------------------------------------------------------------------------------------------
	Column Density calculation
	-----------------------------------------------------------------------------------------------------------------------------------------------
	"""
	itter = 100000
	#define n and nHI and delta r
	n = number_density_H(h100,density, omegab, Xh, ztime)*1e-6 #cm-3
	n_HI = n*H1frac #cm-3
	cm_pos = (posaxis*kPC*1e2)/((1+ztime)*h100)#proper cm
	idx = find_nearest(velaxis, window)
	print(idx)
	delta_r = cm_pos[3] - cm_pos[2]

	#random number seeded so the same random numbers are generated each time
	np.random.seed(50) 
	start_idx = np.random.randint(len(n_HI)-1, size = itter)

		
	
	tile_HI = np.tile(n_HI,2) #incase we sample too close to the end of the array it just cycles round
	tile_pos = np.tile(cm_pos,2*numlos[0])
	tile_Temp = np.tile(temp_1,2)
	tile_density = np.tile(density,2)
	tile_xh1 = np.tile(H1frac,2)
	log_NHI = np.array([])


	logHNI_nothresh = np.array([])
	logHNI_thresh = np.array([])
	av_T = np.array([])
	av_density = np.array([])
	av_xh1 = np.array([])

	for j in range(0,len(start_idx)):
		nHI_slice = tile_HI[start_idx[j]:start_idx[j]+idx]
		xHI_slice = H1frac[start_idx[j]:start_idx[j]+idx]
		pos_slice = tile_pos[start_idx[j]:start_idx[j]+idx]
		av_T = np.append(av_T,np.mean(tile_Temp[start_idx[j]:start_idx[j]+idx]))
		av_density = np.append(av_density,np.mean(tile_density[start_idx[j]:start_idx[j]+idx]))
		av_xh1 = np.append(av_xh1,np.mean(tile_xh1[start_idx[j]:start_idx[j]+idx]))


		#empty array for column desnity
		NHI = integrate.trapezoid(nHI_slice, dx=delta_r) #cm^-2 trap_integ(pos_slice,nHI_slice) #
		logHNI_nothresh = np.append(logHNI_nothresh, np.log10(NHI))

		#N_HI = np.append(N_HI, np.log10(NHI))
	if ii == '10.0':
		d = {'NH1': logHNI_nothresh, 'Temp': av_T, 'rho': av_density, 'xH1': av_xh1}
		data_100 = pd.Series(d)
	if ii == '6.0':
		d = {'NH1': logHNI_nothresh, 'Temp': av_T, 'rho': av_density, 'xH1': av_xh1}
		data_60 = pd.Series(d)
	if ii == '5.4':
		d = {'NH1': logHNI_nothresh, 'Temp': av_T, 'rho': av_density, 'xH1': av_xh1}
		data_54 = pd.Series(d)

#make plot
#Paper specific Matploblib settings
import matplotlib.font_manager
import matplotlib as mpl

column_width=240# call "\the\columnwidth" in LaTeX to find
ppi=72#default ppi, can be left the same

scale=2
fig_width=column_width/ppi*scale#inches
fig_height=3*scale#inches

##SET FONT SIZES
font_small_size = 9
font_medium_size = 12
font_bigger_size = 9

plt.rc('font', size=font_small_size) # controls default text sizes
plt.rc('axes', titlesize=font_small_size) # fontsize of the axes title
plt.rc('axes', labelsize=font_medium_size) # fontsize of the x and y labels
plt.rc('xtick', labelsize=font_small_size) # fontsize of the tick labels
plt.rc('ytick', labelsize=font_small_size) # fontsize of the tick labels
plt.rc('legend', fontsize=font_small_size) # legend fontsize
plt.rc('figure', titlesize=font_bigger_size)


#correct font for MNRAS
#can be found at https://www.fontsquirrel.com/fonts/nimbus-roman-no9-l
#can be installed on Unix systems by putting unzipped folder in directory /home/{user}/.fonts
#run "fc-cache -v" in console to inform system of new font
#print(matplotlib.font_manager.findSystemFonts(fontpaths=None, fontext='ttf'))
#plt.rc('font', family='Nimbus Roman')
# mpl.rcParams["font.family"] = "Nimbus Roman"
# mpl.rcParams['mathtext.fontset'] = 'custom'
# mpl.rcParams['mathtext.rm'] = 'Nimbus Roman'
# mpl.rcParams['mathtext.it'] = 'Nimbus Roman:italic'
# mpl.rcParams['mathtext.bf'] = 'Nimbus Roman:bold'
"""
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
plt.rc('text', usetex=True)
"""

#DPI of MNRAS is 300
mpl.rcParams['figure.dpi'] = 300/scale
# Set up a figure with four panels, with two rows and columns

nrows = 3

ncols = 2

#combinded array for colourbar
combined = np.array([logmean(10**data_100.NH1),  logmean(10**data_60.NH1), logmean(10**data_54.NH1)])
min = np.min(combined)
max = np.max(combined)

from matplotlib.colors import LogNorm
from matplotlib.pyplot import *

#define bin ranges using extent in shape (xmin, xmax, ymin, ymax)
xh1_T = [-5.0,1.0,1.0,5.0]
delta_T = [-1.0,2.0,1.0,5.0]
# axs is a numpy array with dimension (nrows, ncols)
fig, ax = plt.subplots(figsize=(fig_width*2, fig_height*1.5),

                        nrows=nrows,

                        ncols=ncols, sharex='col', sharey ='row')
# Remove horizontal space between axes
levels = 5
fig.subplots_adjust(hspace=0, wspace=0)
ax[0,0].hexbin(np.log10(data_100.rho), np.log10(data_100.Temp), gridsize=50, C= 10**data_100.NH1, reduce_C_function = logmean, vmin=min, vmax=max, cmap='viridis', extent=delta_T)
sc = ax[0,0].hexbin(np.log10(data_100.rho), np.log10(data_100.Temp), gridsize=50, C= 10**data_100.NH1, reduce_C_function = loglen, alpha = 0, extent=delta_T, bins='log')
ax[0,0].set_xlabel(r'$\log_{10} \langle \Delta \rangle$', fontsize = 16)
ax[0,0].set_ylabel(r'$\log_{10} \langle T /K \rangle$', fontsize = 16)
xy = sc.get_offsets()
v = sc.get_array()
print(v)
sns.kdeplot(x=xy[:,0], y=xy[:,1], weights=v, ax=ax[0,0], levels=levels, linestyles='--', linewidths=0.5, color='grey')

sc = ax[0,1].hexbin(np.log10(data_100.xH1), np.log10(data_100.Temp),  gridsize=50, C= 10**data_100.NH1, reduce_C_function =logmean, vmin=min, vmax=max, cmap='viridis', extent=xh1_T) #"""C= 10**data_100.NH1, reduce_C_function = logmean,"""vmin=min, vmax=max, 
ax[0,1].set_xlabel(r'$\log_{10} \langle x_{HI} \rangle$', fontsize = 16)
ax[0,1].text(0.95,0.05, r'$\langle x_{HI} \rangle$ = 0.9',transform=ax[0,0].transAxes, ha='right', va='bottom', fontsize = 12)
sc = ax[0,1].hexbin(np.log10(data_100.xH1), np.log10(data_100.Temp), gridsize=50, C= 10**data_100.NH1, reduce_C_function = loglen, alpha = 0, extent=xh1_T)
xy = sc.get_offsets()
v = sc.get_array()
print(v)
sns.kdeplot(x=xy[:,0], y=xy[:,1], weights=v, ax=ax[0,1], levels=levels, linestyles='--', linewidths=0.5, color='grey')

ax[1,0].hexbin(np.log10(data_60.rho), np.log10(data_60.Temp), gridsize=50, C= 10**data_60.NH1, reduce_C_function = logmean, vmin=min, vmax=max, cmap='viridis', extent=delta_T) #""" C= 10**data_60.NH1, reduce_C_function = logmean,""", vmin=min, vmax=max
ax[1,0].set_xlabel(r'$\log_{10} \langle \Delta \rangle$', fontsize = 16)
ax[1,0].set_ylabel(r'$\log_{10} \langle T /K \rangle$', fontsize = 16)
sc = ax[1,0].hexbin(np.log10(data_60.rho), np.log10(data_60.Temp), gridsize=50, C= 10**data_100.NH1, reduce_C_function = loglen, alpha = 0, extent=delta_T)
xy = sc.get_offsets()
v = sc.get_array()
print(v)
sns.kdeplot(x=xy[:,0], y=xy[:,1], weights=v, ax=ax[1,0], levels=levels, linestyles='--', linewidths=0.5, color='grey')

ax[1,1].hexbin(np.log10(data_60.xH1), np.log10(data_60.Temp), gridsize=50, C= 10**data_60.NH1, reduce_C_function = logmean, cmap='viridis', vmin=min, vmax=max, extent=xh1_T) # """C= 10**data_60.NH1, reduce_C_function = logmean,""", vmin=min, vmax=max
ax[1,1].set_xlabel(r'$\log_{10} \langle x_{HI} \rangle$', fontsize = 16)
ax[1,0].text(0.95,0.05, r'$\langle x_{HI} \rangle = 10^{-1.3}$',transform=ax[1,0].transAxes, ha='right', va='bottom', fontsize=12)
sc = ax[1,1].hexbin(np.log10(data_60.xH1), np.log10(data_60.Temp), gridsize=50, C= 10**data_100.NH1, reduce_C_function = loglen, alpha = 0,  extent=xh1_T) # """C= 10**data_60.NH1, reduce_C_function = logmean,""", vmin=min, vmax=max
xy = sc.get_offsets()
v = sc.get_array()
print(v)
sns.kdeplot(x=xy[:,0], y=xy[:,1], weights=v, ax=ax[1,1], levels=levels, linestyles='--', linewidths=0.5, color='grey')

ax[2,0].hexbin(np.log10(data_54.rho), np.log10(data_54.Temp),gridsize=50, C= 10**data_54.NH1, reduce_C_function = logmean, vmin=min, vmax=max, cmap='viridis', extent=delta_T) # """C= 10**data_54.NH1, reduce_C_function = logmean,""" , vmin=min, vmax=max
ax[2,0].set_xlabel(r'$\log_{10} \langle \Delta \rangle$', fontsize = 16)
ax[2,0].set_ylabel(r'$\log_{10} \langle T /K \rangle$', fontsize = 16)
sc = ax[2,0].hexbin(np.log10(data_54.rho), np.log10(data_54.Temp), gridsize=50, C= 10**data_100.NH1, reduce_C_function = loglen, alpha = 0, extent=delta_T)
xy = sc.get_offsets()
v = sc.get_array()
print(v)
sns.kdeplot(x=xy[:,0], y=xy[:,1], weights=v, ax=ax[2,0], levels=levels, linestyles='--', linewidths=0.5, color='grey')


c = ax[2,1].hexbin(np.log10(data_54.xH1), np.log10(data_54.Temp), gridsize=50, C= 10**data_54.NH1, reduce_C_function = logmean, vmin=min, vmax=max, cmap='viridis', extent=xh1_T) #"""C= 10**data_54.NH1, reduce_C_function = logmean,""" , vmin=min, vmax=max
ax[2,1].set_xlabel(r'$\log_{10} \langle x_{HI} \rangle$', fontsize = 16)
ax[2,0].text(0.95,0.05, r'$\langle x_{HI} \rangle= 10^{-4.2}$',transform=ax[2,0].transAxes, ha='right', va='bottom', fontsize =12)
sc = ax[2,1].hexbin(np.log10(data_54.xH1), np.log10(data_54.Temp), gridsize=50, C= 10**data_100.NH1, reduce_C_function = loglen, alpha = 0, extent=xh1_T) 
xy = sc.get_offsets()
v = sc.get_array()
print(v)
sns.kdeplot(x=xy[:,0], y=xy[:,1], weights=v, ax=ax[2,1], levels=levels, linestyles='--', linewidths=0.5, color='grey')

fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
cb = fig.colorbar(c, cax=cbar_ax)
cb.set_label(label = r'$\langle \log_{10}[N_{HI}/cm^{-2}]\rangle$', size='large')
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.savefig('/home/ppxjf3/hextplot_all_number.pdf')
plt.show()
