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
import random

"""
#argparse variables
parser = argparse.ArgumentParser(description='Files and parameters to calculate mfp.')
parser.add_argument('-f','--folder', help='folder with relevant files', required=True)
parser.add_argument('-z_up','--high_z', help='upper z limit', required=True)
args = parser.parse_args()
"""
#global constants
MPC        = 3.08568025e+22  # m
H0         = 1.0e5/MPC       # s-1
kPC        = const.kpc.value #m
GRAVITY    = const.G.value   #m3 / (kg s2)
PROTONMASS = const.m_p.value #kg
kb = const.k_B.value # J / (K)

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
    return np.log10(np.mean(x))

def jeans_scale(temp, z):
    Y=0.24
    z = float(z)
    
    ne_nH = 1 +  (Y / (4*(1-Y)))
    mu = ((1-Y)*ne_nH+ Y/4)**-1
    T = np.mean(temp)/(10**2)
    lambda_j = 9.8*np.sqrt(10*(1.22/mu)*T*(10/(1+z)))
    return lambda_j


window_prop = 0.5*1000  #proper kpc/h
proper = True
window_func = 'proper'
thresh = 0.5
z = ['6.0', '10.0', '8.0', '7.0', '5.4']
#folder = 'planck1_40_2048_RTzrfit_homog'
#folder = 'planck1_40_2048_ad'
#folder = 'planck1_40_2048_wdm3_RTzrfit'
folder = 'planck1_40_2048_RTzrfit_los'
N_HI_j = np.array([])
N_HI_c = np.array([])
N_HI_p = np.array([])
for jj in z:
	file = '/home/ppxjf3/' + folder + '/los2048_n5000_z' + str(jj)+'00.dat'
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
	--------------------------------------------------------------------------------------------------------------------------------------------
	Column Density calculation
	-----------------------------------------------------------------------------------------------------------------------------------------------
	"""
	itter = 100000
	#define n and nHI and delta r
	n = number_density_H(h100,density, omegab, Xh, ztime)*1e-6 #cm-3
	n_HI = n*H1frac #cm-3
	cm_pos = (posaxis*kPC*1e2)/((1+ztime)*h100)#proper cm
	if window_func == 'proper':
		idx = find_nearest(posaxis/(1+ztime), window_prop)
	if window_func == 'comove':
		window_co = ((1+6.0)/(2))*1000 #cMPc/h
		idx = find_nearest(posaxis, window_co)
	if window_func == 'jeans':
		window_j = jeans_scale(temp_1, ztime)
		idx = find_nearest(posaxis, window_j)
	print(idx)
	delta_r = cm_pos[3] - cm_pos[2]

		#random number seeded so the same random numbers are generated each time
	np.random.seed(50) 
	start_idx = np.random.randint(len(n_HI)-1, size = itter)

			
		
	tile_HI = np.tile(n_HI,2) #incase we sample too close to the end of the array it just cycles round
	tile_pos = np.tile(cm_pos,2*numlos[0])
	tile_Temp = np.tile(temp_1,2)
	tile_density = np.tile(density,2)
	log_NHI = np.array([])


	logHNI_nothresh = np.array([])
	logHNI_thresh = np.array([])
	av_T = np.array([])
	av_density = np.array([])
	log_av_T = np.array([])
	log_av_density = np.array([])

	for j in range(0,len(start_idx)):
		nHI_slice = tile_HI[start_idx[j]:start_idx[j]+idx]
		xHI_slice = H1frac[start_idx[j]:start_idx[j]+idx]
		pos_slice = tile_pos[start_idx[j]:start_idx[j]+idx]
		av_T = np.append(av_T,np.mean(tile_Temp[start_idx[j]:start_idx[j]+idx]))
		av_density = np.append(av_density,np.mean(tile_density[start_idx[j]:start_idx[j]+idx]))
		log_av_T = np.append(log_av_T,np.mean(np.log10(tile_Temp[start_idx[j]:start_idx[j]+idx])))
		log_av_density = np.append(log_av_density,np.mean(np.log10(tile_density[start_idx[j]:start_idx[j]+idx])))


		#empty array for column desnity
		NHI = integrate.trapezoid(nHI_slice, dx=delta_r) #cm^-2 trap_integ(pos_slice,nHI_slice) #
		logHNI_nothresh = np.append(logHNI_nothresh, np.log10(NHI))
		if np.sum(xHI_slice) < thresh:
			
			logHNI_thresh = np.append(logHNI_thresh, np.log10(NHI))
		else:
			continue
		#N_HI = np.append(N_HI, np.log10(NHI))
	if jj == '10.0':
		NHI_100 = logHNI_nothresh
		NHI_100_thresh = logHNI_thresh
		av_xHI_100 = np.mean(H1frac)
		dr_100 = normalisation_factor(velaxis[idx], itter, h100) #MPc
	if jj == '8.0':
		NHI_80 = logHNI_nothresh
		NHI_80_thresh = logHNI_thresh
		av_xHI_80 = np.mean(H1frac)
		dr_80 = normalisation_factor(velaxis[idx], itter, h100) #MPc
	if jj == '7.0':
		NHI_70 = logHNI_nothresh
		NHI_70_thresh = logHNI_thresh
		av_xHI_70 = np.mean(H1frac)
		dr_70 = normalisation_factor(velaxis[idx], itter, h100) #MPc
	if jj == '6.0':
		NHI_60 = logHNI_nothresh
		NHI_60_thresh = logHNI_thresh
		av_xHI_60 = np.mean(H1frac)
		dr_60 = normalisation_factor(velaxis[idx], itter, h100) #MPc
	if jj == '5.4':
		NHI_54 = logHNI_nothresh
		NHI_54_thresh = logHNI_thresh
		av_xHI_54 = np.mean(H1frac)
		dr_54 = normalisation_factor(velaxis[idx], itter, h100) #MPc


pdf_100, bins_100, centers_100 = pdf_calc(np.min(NHI_100), np.max(NHI_100), 18, NHI_100, dr_100)
pdf_80, bins_80, centers_80 = pdf_calc(np.min(NHI_80), np.max(NHI_80), 18, NHI_80, dr_80)
pdf_70, bins_70, centers_70 = pdf_calc(np.min(NHI_70), np.max(NHI_70), 18, NHI_70, dr_70)
pdf_60, bins_60, centers_60 = pdf_calc(np.min(NHI_60), np.max(NHI_60), 18, NHI_60, dr_60)
pdf_54, bins_54, centers_54 = pdf_calc(np.min(NHI_54), np.max(NHI_54), 18, NHI_54, dr_54)

pdf_100_thresh, bins_100_thresh, centers_100_thresh = pdf_calc(np.min(NHI_100_thresh), np.max(NHI_100_thresh), 18, NHI_100_thresh, dr_100)
pdf_80_thresh, bins_80_thresh, centers_80_thresh = pdf_calc(np.min(NHI_80_thresh), np.max(NHI_80_thresh), 18, NHI_80_thresh, dr_80)
pdf_70_thresh, bins_70_thresh, centers_70_thresh = pdf_calc(np.min(NHI_70_thresh), np.max(NHI_70_thresh), 18, NHI_70_thresh, dr_70)
pdf_60_thresh, bins_60_thresh, centers_60_thresh = pdf_calc(np.min(NHI_60_thresh), np.max(NHI_60_thresh), 18, NHI_60_thresh, dr_60)
pdf_54_thresh, bins_54_thresh, centers_54_thresh = pdf_calc(np.min(NHI_54_thresh), np.max(NHI_54_thresh), 18, NHI_54_thresh, dr_54)
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

nrows = 1

ncols = 1

# axs is a numpy array with dimension (nrows, ncols)
colour = ['b', 'm', 'darkorange', 'g', 'r']
fig, ax = plt.subplots(figsize=(fig_width, fig_height),

                        nrows=nrows,

                        ncols=ncols)


ax.plot(centers_100, np.log10(pdf_100), color = colour[0], label=r'$\langle x_{\mathrm{H1}} \rangle = 0.9$')
ax.plot(centers_80, np.log10(pdf_80), color = colour[1], label=r'$\langle x_{\mathrm{H1}} \rangle = 10^{-0.2}$')
ax.plot(centers_70, np.log10(pdf_70), color = colour[2], label=r'$\langle x_{\mathrm{H1}} \rangle = 10^{-0.4}$')
ax.plot(centers_60, np.log10(pdf_60), color = colour[3], label=r'$\langle x_{\mathrm{H1}} \rangle = 10^{-1.3}$')
ax.plot(centers_54, np.log10(pdf_54), color = colour[4], label=r'$\langle x_{\mathrm{H1}} \rangle = 10^{-4.2}$')

ax.plot(centers_54_thresh, np.log10(pdf_54_thresh), linestyle = '--', color = colour[4])#, label=r'$\langle x_{\mathrm{H1}} \rangle = 10^{-4.2}$')#.format(np.log10(av_xHI_54)))
ax.plot(centers_60_thresh, np.log10(pdf_60_thresh), linestyle='--', color = colour[3])#, label=r'$\langle x_{\mathrm{H1}} \rangle = 10^{-1.3}$')#label='xHI=1e{:.2f}'.format(np.log10(av_xHI_60)))
ax.plot(centers_70_thresh, np.log10(pdf_70_thresh), color = colour[2], linestyle='--')#, label=r'$\langle x_{\mathrm{H1}} \rangle = 10^{-0.4}$')#label='xHI=1e{:.2f}'.format(np.log10(av_xHI_70)))
ax.plot(centers_80_thresh, np.log10(pdf_80_thresh),linestyle = '--', color = colour[1])#, label=r'$\langle x_{\mathrm{H1}} \rangle = 10^{-0.2}$')#label='xHI=1e{:.2f}'.format(np.log10(av_xHI_80)))
ax.plot(centers_100_thresh, np.log10(pdf_100_thresh), color = colour[0], linestyle = '--')

ax.set_xlabel(r'$log_{10}[N_{\mathrm{H1}} / \mathrm{cm}^{-2}]$', fontsize = 18)
ax.set_ylabel(r'$log_{10}[\frac{ \partial ^2n}{ \partial log_{10}[N_{\mathrm{H1}}] \partial R} /\mathrm{ pMPc}^{-1}]$', fontsize = 18)
#ax[0].set_xticks(fontsize=16)
#ax[0].set_yticks(fontsize=16)
ax.legend(frameon=False, fontsize = 11)
plt.tight_layout()

base = '/home/ppxjf3/planck1_40_2048_RTzrfit_los/column_desnity_data/'

plt.savefig(base + 'paper1_fig3_window_fixed' +str(window_func)+'.pdf', bbox_inches='tight', pad_inches=0.02)  # remove whitespace
plt.show()