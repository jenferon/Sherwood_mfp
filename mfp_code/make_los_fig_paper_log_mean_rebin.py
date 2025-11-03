import numpy as np
import matplotlib.pyplot as plt
import math
from scipy import interpolate
#import seaborn as sns
from mpl_toolkits.axes_grid1 import make_axes_locatable
from astropy import constants as const
import random
import scipy
import os
import pandas as pd
#Paper specific Matploblib settings
import matplotlib.font_manager
import matplotlib as mpl
import seaborn as sns
from scipy.stats import bootstrap
from scipy import stats
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

import matplotlib
matplotlib.style.use('/home/ppxjf3/paper_params.mplstyle')

#global constants
MPC        = 3.08568025e+22  # m
H0         = 1.0e5/MPC       # s-1
kPC        = const.kpc.value #m
GRAVITY    = const.G.value   #m3 / (kg s2)
PROTONMASS = const.m_p.value #kg
M_sol = const.M_sun.value #kg
cross_sec = 6.3*10**-18 #cm2 for 912 amstrongs
itter = 100000



class data_file(object):
	def __init__(self, filename, filename_halo):
		#print(filename)
		readdata = open(filename, 'rb')
		# Header data
		self.ztime  = np.fromfile(readdata,dtype=np.double,count=1) # redshift
		self.omegaM = np.fromfile(readdata,dtype=np.double,count=1) # Omega_m (matter density)
		self.omegaL = np.fromfile(readdata,dtype=np.double,count=1) # Omega_L (Lambda density)
		self.omegab = np.fromfile(readdata,dtype=np.double,count=1) # Omega_b (baryon density)
		self.h100   = np.fromfile(readdata,dtype=np.double,count=1) # Hubble constant, H0 / 100 km/s/Mpc
		self.box100 = np.fromfile(readdata,dtype=np.double,count=1) # Box size in comoving kpc/h
		self.Xh     = np.fromfile(readdata,dtype=np.double,count=1) # Hydrogen fraction by mass
		self.nbins  = np.fromfile(readdata,dtype=np.int32,count=1)  # Number of pixels in each line of sight
		self.numlos = np.fromfile(readdata,dtype=np.int32,count=1)  # Number of lines of sight
		print(self.numlos)

		# Line of sight locations in box 
		self.iaxis  = np.fromfile(readdata,dtype=np.int32,count=self.numlos[0])  # projection axis, x=1, y=2, z=3
		self.xaxis  = np.fromfile(readdata,dtype=np.double,count=self.numlos[0]) # x-coordinate in comoving kpc/h
		self.yaxis  = np.fromfile(readdata,dtype=np.double,count=self.numlos[0]) # y-coordinate in comoving kpc/h
		self.zaxis  = np.fromfile(readdata,dtype=np.double,count=self.numlos[0]) # z-coordinate in comoving kpc/h

		# Line of sight scale
		self.posaxis = np.fromfile(readdata,dtype=np.double,count=self.nbins[0]) # comoving kpc/h
		self.posaxis_p = self.posaxis/((1+self.ztime)*self.h100*1000) #pMpc
		self.velaxis = np.fromfile(readdata,dtype=np.double,count=self.nbins[0]) # km/s hubble velocity
		
		# Gas density, rho/<rho>
		self.density = np.fromfile(readdata,dtype=np.double,count=self.nbins[0]*self.numlos[0])

		# H1 fraction, fH1 = nH1/nH
		self.H1frac  = np.fromfile(readdata,dtype=np.double,count=self.nbins[0]*self.numlos[0])

		# Temperature, K
		self.temp    = np.fromfile(readdata,dtype=np.double,count=self.nbins[0]*self.numlos[0])

		# Peculiar velocity, km/s
		self.vpec    = np.fromfile(readdata,dtype=np.double,count=self.nbins[0]*self.numlos[0])

		self.Gammaker = np.fromfile(readdata,dtype=np.double,count=self.nbins[0]*self.numlos[0]) #/* Gamma_HIc [s^-1] */
		# Close the binary file
		readdata.close()

		nhaloes = 3000

        # Read in halo list with CoM and mass
		with open(filename_halo, 'r') as file:
			lines = file.readlines()

		iaxis_fof = np.zeros(nhaloes, dtype=int)  # axis over which we want to draw the FoF halo
		group_x = np.zeros(nhaloes, dtype=float)  # x CoM ckpc/h
		group_y = np.zeros(nhaloes, dtype=float)  # y CoM ckpc/h
		group_z = np.zeros(nhaloes, dtype=float)  # z CoM ckpc/h
		mass_tot = np.zeros(nhaloes, dtype=float)  # log10(M / Msol/h)

        # Read in data
		count = 0
		for line in lines:
			in1, in2, in3, in4, in5 = map(float, line.strip().split())
			iaxis_fof[count] = int(in1)
			group_x[count] = in2
			group_y[count] = in3
			group_z[count] = in4
			mass_tot[count] = in5
			count += 1

		self.i_axis = iaxis_fof
		self.group_x = group_x
		self.group_y = group_y
		self.group_z = group_z
		self.mass_tot = mass_tot

	def number_density_H(self):
	#Get nH(z) at the critical density in proper cgs units

		rhoc = 3.0 * (H0*self.h100)**2/(8.0 * np.pi * GRAVITY) #kg m-3

		self.n_H = rhoc*self.omegab*(self.density)*self.Xh*(1.0 + self.ztime)**3.0 / PROTONMASS #m-3
  
		self.n_HI = self.n_H*self.H1frac #m-3
  
	def unit_conversion(self):
		cm_pos = (self.posaxis*kPC*1e2)/((1+self.ztime)*self.h100)#proper cm
		delta_r = abs(cm_pos[300] - cm_pos[299])#
		return cm_pos, delta_r

	def los_calc(self,k):
		los = self.density[k*self.nbins[0] : (k+1)*self.nbins[0]]
		los_n_HI = self.n_HI[k*self.nbins[0] : (k+1)*self.nbins[0]]
		los_x_HI = self.H1frac[k*self.nbins[0] : (k+1)*self.nbins[0]]
		los_n_H = self.n_H[k*self.nbins[0] : (k+1)*self.nbins[0]]
		return los, los_n_HI, los_x_HI, los_n_H

	def find_nearest(self, value):
		array = np.asarray(self.posaxis)
		idx = (np.abs(self.posaxis - value)).argmin()
		return idx

	def r_vir_calc(self, M):
		#print(M)
		omegam_z = (self.omegaM * (1.0 + self.ztime)**3.0) / (self.omegaM * (1.0 + self.ztime)**3.0 + self.omegaL)

		d = omegam_z - 1.0

		Delta_c = 18.0 * np.pi**2

		Delta_cz = Delta_c + 82.0 * d - 39.0 * d**2

		rvir = 0.784 * np.cbrt((10**M) / 1.0e8)* (((self.omegaM * Delta_cz) / (omegam_z * Delta_c))**(-1.0 / 3.0)) * (10.0 / (1.0 + self.ztime))  # pkpc/h

		"""omegaK = 1 - self.omegaM -self.omegaL - self.omegab
		omegaM_z = (self.omegaM*(1+self.ztime)**3)/(self.omegaM*(1+self.ztime)**3 + self.omegaL)# + omegaK*(1+self.ztime)**2)
		d = omegaM_z - 1
		delta_c = 18*np.pi**2 + 82*d - 39*d**2
		r_vir_proper = 0.784*np.cbrt(M/(1.0e8))*np.cbrt((omegaM_z*18*np.pi**2)/(self.omegaM*delta_c))*(10/(1+self.ztime)) #proper kpc/h"""
		r_vir = rvir/(self.h100) #proper kpc
		#print(self.posaxis[0]-self.posaxis[1])
		steps = r_vir/(self.posaxis[0]-self.posaxis[1])
		steps = np.abs(round(steps[0]))
		#print('virial radius =' +str(r_vir)+ ' (c kpc/h) for mass '+str(M))
		#print(steps)
		
		return r_vir*10**-3 #pMpc
    
	def get_los(self):
        
		self.number_density_H()

		#empty arrays
		los_nH = np.empty([int(self.numlos*2),self.nbins[0]])
		los_nHI = np.empty([int(self.numlos*2),self.nbins[0]])
		los_xHI = np.empty([int(self.numlos*2),self.nbins[0]])
		delta = np.empty([int(self.numlos*2),self.nbins[0]])
		r_vir = np.zeros([int(self.numlos)-1])
		ii = 0
		jj=0

		for k in range(0, int(self.numlos)-1):

			los_rho, los_n_HI, los_x_HI, los_n_H = self.los_calc(k)
	
			if self.iaxis[k] == 1:
				max_ind = self.find_nearest(self.group_x[k])
				r_vir[ii] = self.r_vir_calc(self.mass_tot[k])

			if self.iaxis[k] == 2:
				max_ind = self.find_nearest(self.group_y[k])
				r_vir[ii] = self.r_vir_calc(self.mass_tot[k])
			
			if self.iaxis[k] == 3:
				max_ind = self.find_nearest(self.group_z[k])
				r_vir[ii] = self.r_vir_calc(self.mass_tot[k])

			
			cut_los = np.roll(los_rho, -max_ind)
			cut_n_HI = np.roll(los_n_HI, -max_ind)
			cut_nH =  np.roll(los_n_H, -max_ind)
			cut_XHI =  np.roll(los_x_HI, -max_ind)
			
			los_nH[jj,:] = cut_nH
			los_nHI[jj,:] = cut_n_HI
			los_xHI[jj,:] = cut_XHI
			delta[jj,:] = cut_los
			jj += 1
			los_nH[jj,:] = cut_nH[::-1]
			los_nHI[jj,:] = cut_n_HI[::-1]
			los_xHI[jj,:] = cut_XHI[::-1]
			delta[jj,:] = cut_los[::-1]
			jj += 1
			ii += 1

		self.delta = delta
		self.xHI = los_xHI
		self.nHI = los_nHI
		#print(r_vir_idx)
		self.rvir = np.mean(r_vir)
		print(self.rvir)
 
	def get_los_rand(self):
            
		self.number_density_H()
            
		np.random.seed(767) 
		start_value = np.random.randint(self.numlos-1, size = itter)
		print(self.ztime)
		print(np.mean(self.H1frac))

		#empty arrays
		los_nH = np.empty([itter,self.nbins[0]])
		los_nHI = np.empty([itter,self.nbins[0]])
		los_xHI = np.empty([itter,self.nbins[0]])
		delta = np.empty([itter,self.nbins[0]])
		o = 0
        
		for k in start_value:

			los_rho, los_n_HI, los_x_HI, los_n_H = self.los_calc(k)
                
			los_nH[o,:] = (los_n_H)
			los_nHI[o,:] = (los_n_HI)
			los_xHI[o,:] = (los_x_HI)
			delta[o,:] = (los_rho)
			o += 1
		self.delta = delta
		self.xHI = los_xHI
		self.nHI = los_nHI

	def rebin_new(self, avg_density, avg_H1frac, scale=4):
		nbins_new = int(len(self.posaxis)/scale)
		print(nbins_new)
		dR            = self.box100/nbins_new
		halfdR        = dR/2
		posaxis_rebin = min(self.posaxis) + dR*np.arange(nbins_new)
		avg_density_rebin = np.empty([int(self.numlos*2),int(nbins_new)])
		avg_H1frac_rebin = np.empty([int(self.numlos*2),int(nbins_new)])
		for i in range(0, int(self.numlos*2-1)):
			for j in range(0, int(nbins_new-1)):
				bin_lower = posaxis_rebin[j] - halfdR
				bin_upper = posaxis_rebin[j] + halfdR
				indbin = np.argwhere((bin_lower < self.posaxis) & (self.posaxis < bin_upper))
				temp_HI_frac = np.empty([int(len(self.posaxis)),len(indbin)])
				temp_density = np.empty([int(len(self.posaxis)),len(indbin)])
				print(indbin)
				print(i)
				print(self.numlos)
				for k in range(0,len(indbin)):
					temp_density[:,k] = avg_density[i,indbin[k]]
					temp_HI_frac[:,k] = avg_H1frac[i,indbin[k]]
				avg_density_rebin[i,j] = np.mean(temp_density)
				avg_H1frac_rebin[i,j] = np.mean(temp_HI_frac)

		RpMpc_rebin = 10**3*posaxis_rebin/(self.h100*(1.0+self.ztime))
		return RpMpc_rebin, np.mean(avg_density_rebin, axis=0), np.mean(avg_H1frac_rebin, axis=0)

	def rebin(self, scale=2):
		nbins_new = len(self.posaxis)//scale

		avg_H1frac = np.mean(self.xHI, axis=0)
		avg_density = np.mean(self.delta, axis=0)

		dR            = self.box100/nbins_new
		halfdR        = dR/2
		RpMpc_rebin = min(self.posaxis) + dR*np.arange(nbins_new)
		avg_density_rebin = np.zeros([nbins_new])
		avg_H1frac_rebin = np.zeros([nbins_new])
		for j in range(0, int(nbins_new)):
			bin_lower = RpMpc_rebin[j] - halfdR
			bin_upper = RpMpc_rebin[j] + halfdR
			indbin = np.where((self.posaxis >= bin_lower) & (self.posaxis < bin_upper))
			print(indbin)
			if len(indbin) == 0:
				print('!! Error -- too few pixels for rebinning !!')
				raise ValueError('Too few pixels for rebinning')
			else:
				avg_density_rebin[j] = np.mean(avg_density[indbin])
				avg_H1frac_rebin[j] = np.mean(avg_H1frac[indbin])
		RpMpc_rebin =  1.0e-3 * RpMpc_rebin / (self.h100 * (1.0 + self.ztime))
		return RpMpc_rebin, avg_density_rebin, avg_H1frac_rebin
base = '/home/ppxjf3/planck1_40_2048_RTzrfit_haloes/'

#read in data
mfp_class10_1 = data_file(str(base + 'halolow_los2048_n3000_z10.000.dat'), str(base + 'haloes_low3000_z10.000.dat'))
mfp_class10_1.get_los()
RpMpc_rebin10_1, avg_density_rebin10_1, avg_H1frac_rebin10_1 = mfp_class10_1.rebin()
mfp_class10_2 = data_file(str(base + 'halomid_los2048_n3000_z10.000.dat'), str(base + 'haloes_mid3000_z10.000.dat'))
mfp_class10_2.get_los()
RpMpc_rebin10_2, avg_density_rebin10_2, avg_H1frac_rebin10_2 = mfp_class10_2.rebin()
mfp_class10_3 = data_file(str(base + 'halohigh_los2048_n3000_z10.000.dat'), str(base + 'haloes_high3000_z10.000.dat'))
mfp_class10_3.get_los()
RpMpc_rebin10_3, avg_density_rebin10_3, avg_H1frac_rebin10_3 = mfp_class10_3.rebin()
mfp_class10_rand = data_file('/home/ppxjf3/planck1_40_2048_RTzrfit_los/los2048_n5000_z10.000.dat', str(base + 'haloes_high3000_z10.000.dat'))
mfp_class10_rand.get_los_rand()

mfp_class6_1 = data_file(str(base + 'halolow_los2048_n3000_z6.000.dat'), str(base + 'haloes_low3000_z6.000.dat'))
mfp_class6_1.get_los()
RpMpc_rebin6_1, avg_density_rebin6_1, avg_H1frac_rebin6_1 = mfp_class6_1.rebin()
mfp_class6_2 = data_file(str(base + 'halomid_los2048_n3000_z6.000.dat'), str(base + 'haloes_mid3000_z6.000.dat'))
mfp_class6_2.get_los()
RpMpc_rebin6_2, avg_density_rebin6_2, avg_H1frac_rebin6_2 = mfp_class6_2.rebin()
mfp_class6_3 = data_file(str(base +'halohigh_los2048_n3000_z6.000.dat'), str(base + 'haloes_high3000_z6.000.dat'))
mfp_class6_3.get_los()
RpMpc_rebin6_3, avg_density_rebin6_3, avg_H1frac_rebin6_3 = mfp_class6_3.rebin()
mfp_class6_rand = data_file('/home/ppxjf3/planck1_40_2048_RTzrfit_los/los2048_n5000_z6.000.dat', str(base + 'haloes_high3000_z6.000.dat'))
mfp_class6_rand.get_los_rand()

mfp_class54_1 = data_file(str(base + 'halolow_los2048_n3000_z5.400.dat'), str(base + 'haloes_low3000_z5.400.dat'))
mfp_class54_1.get_los()
RpMpc_rebin54_1, avg_density_rebin54_1, avg_H1frac_rebin54_1 = mfp_class54_1.rebin()
mfp_class54_2 = data_file(str(base + 'halomid_los2048_n3000_z5.400.dat'), str(base + 'haloes_mid3000_z5.400.dat'))
mfp_class54_2.get_los()
RpMpc_rebin54_2, avg_density_rebin54_2, avg_H1frac_rebin54_2 = mfp_class54_2.rebin()
mfp_class54_3 = data_file(str(base +'halohigh_los2048_n3000_z5.400.dat'), str(base + 'haloes_high3000_z5.400.dat'))
mfp_class54_3.get_los()
RpMpc_rebin54_3, avg_density_rebin54_3, avg_H1frac_rebin54_3 = mfp_class54_3.rebin()
mfp_class54_rand = data_file('/home/ppxjf3/planck1_40_2048_RTzrfit_los/los2048_n5000_z5.400.dat', str(base + 'haloes_high3000_z5.400.dat'))
mfp_class54_rand.get_los_rand()


column_width=240# call "\the\columnwidth" in LaTeX to find
ppi=72#default ppi, can be left the same

#plot graph
scale=2
fig_width=column_width/ppi*scale#inches
fig_height=3*scale#inches


##SET FONT SIZES
font_small_size = 9
font_medium_size = 12
font_bigger_size = 14

plt.rc('font', size=font_small_size) # controls default text sizes
plt.rc('axes', titlesize=font_small_size) # fontsize of the axes title
plt.rc('axes', labelsize=font_medium_size) # fontsize of the x and y labels
plt.rc('xtick', labelsize=font_bigger_size) # fontsize of the tick labels
plt.rc('ytick', labelsize=font_bigger_size) # fontsize of the tick labels
plt.rc('legend', fontsize=font_small_size) # legend fontsize
plt.rc('figure', titlesize=font_bigger_size)

#DPI of MNRAS is 300
mpl.rcParams['figure.dpi'] = 300/scale
# Set up a figure with four panels, with two rows and columns

nrows = 3

ncols = 2

#make graph
ymax=0.25
fig, ax = plt.subplots(figsize=(fig_width * 2, fig_height * 1.5),

                        nrows=nrows,

                        ncols=ncols, sharex= 'col')
ax[0,0].plot(RpMpc_rebin10_3, np.log10(avg_density_rebin10_3), c='cyan', label=r'$\langle \rm M/M_\odot h^{-1} \rangle=10^{9.9}$')
ax[0,0].axvline(mfp_class10_3.rvir, c='cyan', linestyle=':', ymin=0 ,ymax=ymax, lw=2)
ax[0,0].axvline(10*mfp_class10_3.rvir, c='cyan', linestyle='dashed', ymin=0 ,ymax=ymax, lw=2)
ax[0,0].plot(RpMpc_rebin10_2, np.log10(avg_density_rebin10_2), c='purple', label=r'$\langle \rm M/M_\odot h^{-1} \rangle=10^{9.4}$')
ax[0,0].axvline(mfp_class10_2.rvir, c='purple', linestyle=':', ymin=0 ,ymax=ymax, lw=2)
ax[0,0].axvline(10*mfp_class10_2.rvir, c='purple', linestyle='dashed', ymin=0 ,ymax=ymax, lw=2)
ax[0,0].plot(RpMpc_rebin10_1, np.log10(avg_density_rebin10_1), c='red', label=r'$\langle \rm M/M_\odot h^{-1} \rangle=10^{9.1}$')
ax[0,0].axvline(mfp_class10_1.rvir, c='red', linestyle=':', ymin=0 ,ymax=ymax, lw=2)
ax[0,0].axvline(10*mfp_class10_1.rvir, c='red', linestyle='dashed', ymin=0 ,ymax=ymax, lw=2)
ax[0,0].plot(mfp_class10_rand.posaxis_p, np.log10(np.mean(mfp_class10_rand.delta, axis=0)), c='green')
ax[0,0].set_ylabel(r'$ \log_{10}\langle \Delta \rangle $',fontsize=16) 
ax[0,0].legend(frameon=False, fontsize = 14)

ax[1,0].plot(RpMpc_rebin6_3, np.log10(avg_density_rebin6_3), c='cyan', label=r'$\langle \rm M/M_\odot h^{-1} \rangle=10^{10.7}$')
ax[1,0].plot(RpMpc_rebin6_2, np.log10(avg_density_rebin6_2), c='purple', label=r'$\langle \rm M/M_\odot h^{-1} \rangle=10^{10.0}$')
ax[1,0].plot(RpMpc_rebin6_1, np.log10(avg_density_rebin6_1), c='red', label=r'$\langle \rm M/M_\odot h^{-1} \rangle=10^{9.3}$')
ax[1,0].plot(mfp_class6_rand.posaxis_p, np.log10(np.mean(mfp_class6_rand.delta, axis=0)), c='green')
ax[1,0].axvline(mfp_class6_3.rvir, c='cyan', linestyle=':', ymin=0 ,ymax=ymax, lw=2)
ax[1,0].axvline(mfp_class6_2.rvir, c='purple', linestyle=':', ymin=0 ,ymax=ymax, lw=2)
ax[1,0].axvline(mfp_class6_1.rvir, c='red', linestyle=':', ymin=0 ,ymax=ymax, lw=2)
ax[1,0].axvline(10*mfp_class6_3.rvir, c='cyan', linestyle='dashed', ymin=0 ,ymax=ymax, lw=2)
ax[1,0].axvline(10*mfp_class6_2.rvir, c='purple', linestyle='dashed', ymin=0 ,ymax=ymax, lw=2)
ax[1,0].axvline(10*mfp_class6_1.rvir, c='red', linestyle='dashed', ymin=0 ,ymax=ymax, lw=2)
ax[1,0].set_ylabel(r'$ \log_{10} \langle\Delta \rangle $',fontsize=16)
ax[1,0].legend(frameon=False, fontsize = 14)

ax[2,0].plot(RpMpc_rebin54_3, np.log10(avg_density_rebin54_3), c='cyan', label=r'$\langle \rm M/M_\odot h^{-1} \rangle=10^{10.9}$')
ax[2,0].plot(RpMpc_rebin54_2, np.log10(avg_density_rebin54_2), c='purple', label=r'$\langle \rm M/M_\odot h^{-1} \rangle=10^{10.1}$')
ax[2,0].plot(RpMpc_rebin54_1, np.log10(avg_density_rebin54_1), c='red', label=r'$\langle \rm M/M_\odot h^{-1} \rangle=10^{9.3}$')
ax[2,0].plot(mfp_class54_rand.posaxis_p, np.log10(np.mean(mfp_class54_rand.delta, axis=0)), c='green')
ax[2,0].axvline(mfp_class54_3.rvir, c='cyan', linestyle=':', ymin=0 ,ymax=ymax, lw=2)
ax[2,0].axvline(mfp_class54_2.rvir, c='purple', linestyle=':', ymin=0 ,ymax=ymax, lw=2)
ax[2,0].axvline(mfp_class54_1.rvir, c='red', linestyle=':', ymin=0 ,ymax=ymax, lw=2)
ax[2,0].axvline(10*mfp_class54_3.rvir, c='cyan', linestyle='dashed', ymin=0 ,ymax=ymax, lw=2)
ax[2,0].axvline(10*mfp_class54_2.rvir, c='purple', linestyle='dashed', ymin=0 ,ymax=ymax, lw=2)
ax[2,0].axvline(10*mfp_class54_1.rvir, c='red', linestyle='dashed', ymin=0 ,ymax=ymax, lw=2)
ax[2,0].set_ylabel(r'$\log_{10} \langle \Delta \rangle $',fontsize=16)
ax[2,0].set_xlabel('R [pMpc]',fontsize=16)
ax[2,0].legend(frameon=False, fontsize = 14)



ax[0,1].plot(RpMpc_rebin10_3, np.log10(avg_H1frac_rebin10_3), c='cyan', label=r'$\langle \rm M/M_\odot h^{-1} \rangle=10^{9.8}$')
ax[0,1].plot(RpMpc_rebin10_2, np.log10(avg_H1frac_rebin10_2), c='purple', label=r'$\langle \rm M/M_\odot h^{-1} \rangle=10^{9.4}$')
ax[0,1].plot(RpMpc_rebin10_1, np.log10(avg_H1frac_rebin10_1), c='red', label=r'$\langle \rm M/M_\odot h^{-1} \rangle=10^{9.1}$')
ax[0,1].plot(mfp_class10_rand.posaxis_p, np.log10(np.mean(mfp_class10_rand.xHI, axis=0)), c='green')
ax[0,1].axvline(mfp_class10_3.rvir, c='cyan', linestyle=':', ymin=0 ,ymax=ymax, lw=2)
ax[0,1].axvline(mfp_class10_2.rvir, c='purple', linestyle=':', ymin=0 ,ymax=ymax, lw=2)
ax[0,1].axvline(mfp_class10_1.rvir, c='red', linestyle=':', ymin=0 ,ymax=ymax, lw=2)
ax[0,1].axvline(10*mfp_class10_3.rvir, c='cyan', linestyle='dashed', ymin=0 ,ymax=ymax, lw=2)
ax[0,1].axvline(10*mfp_class10_2.rvir, c='purple', linestyle='dashed', ymin=0 ,ymax=ymax, lw=2)
ax[0,1].axvline(10*mfp_class10_1.rvir, c='red', linestyle='dashed', ymin=0 ,ymax=ymax, lw=2)
ax[0,1].set_ylabel(r'$ \log_{10}\langle x_{\rm HI} \rangle $',fontsize=16) 
ax[0,1].text(0.95,0.15, r'$\langle x_{\rm HI} \rangle = 0.90$  $(z=10.0)$',transform=ax[0,1].transAxes, ha='right', va='bottom', fontsize=14)

ax[1,1].plot(RpMpc_rebin6_3, np.log10(avg_H1frac_rebin6_3), c='cyan')
ax[1,1].plot(RpMpc_rebin6_2, np.log10(avg_H1frac_rebin6_2), c='purple')
ax[1,1].plot(RpMpc_rebin6_1, np.log10(avg_H1frac_rebin6_1), c='red')
ax[1,1].plot(mfp_class6_rand.posaxis_p, np.log10(np.mean(mfp_class6_rand.xHI, axis=0)), c='green')
ax[1,1].axvline(mfp_class6_3.rvir, c='cyan', linestyle=':', ymin=0 ,ymax=ymax, lw=2)
ax[1,1].axvline(mfp_class6_2.rvir, c='purple', linestyle=':', ymin=0 ,ymax=ymax, lw=2)
ax[1,1].axvline(mfp_class6_1.rvir, c='red', linestyle=':', ymin=0 ,ymax=ymax, lw=2)
ax[1,1].axvline(10*mfp_class6_3.rvir, c='cyan', linestyle='dashed', ymin=0 ,ymax=ymax, lw=2)
ax[1,1].axvline(10*mfp_class6_2.rvir, c='purple', linestyle='dashed', ymin=0 ,ymax=ymax, lw=2)
ax[1,1].axvline(10*mfp_class6_1.rvir, c='red', linestyle='dashed', ymin=0 ,ymax=ymax, lw=2)
ax[1,1].set_ylabel(r'$\log_{10}\langle  x_{\rm HI} \rangle $',fontsize=16)
ax[1,1].text(0.95,0.5, r'$\langle x_{\rm HI} \rangle = 0.05$  $(z=6.0)$',transform=ax[1,1].transAxes, ha='right', va='bottom', fontsize=14)

ax[2,1].plot(RpMpc_rebin54_3, np.log10(avg_H1frac_rebin54_3), c='cyan')
ax[2,1].plot(RpMpc_rebin54_2, np.log10(avg_H1frac_rebin54_2), c='purple')
ax[2,1].plot(RpMpc_rebin54_1, np.log10(avg_H1frac_rebin54_1), c='red')
ax[2,1].plot(mfp_class54_rand.posaxis_p, np.log10(np.mean(mfp_class54_rand.xHI, axis=0)), c='green')
ax[2,1].axvline(mfp_class54_3.rvir, c='cyan', linestyle=':', ymin=0 ,ymax=ymax, lw=2)
ax[2,1].axvline(mfp_class54_2.rvir, c='purple', linestyle=':', ymin=0 ,ymax=ymax, lw=2)
ax[2,1].axvline(mfp_class54_1.rvir, c='red', linestyle=':', ymin=0 ,ymax=ymax, lw=2)
ax[2,1].axvline(10*mfp_class54_3.rvir, c='cyan', linestyle='dashed', ymin=0 ,ymax=ymax, lw=2)
ax[2,1].axvline(10*mfp_class54_2.rvir, c='purple', linestyle='dashed', ymin=0 ,ymax=ymax, lw=2)
ax[2,1].axvline(10*mfp_class54_1.rvir, c='red', linestyle='dashed', ymin=0 ,ymax=ymax, lw=2)
ax[2,1].set_ylabel(r'$\log_{10}\langle  x_{\rm HI} \rangle $',fontsize=16)
ax[2,1].set_xlabel('R [pMpc]',fontsize=16)
ax[2,1].text(0.95,0.45, r'$\langle x_{\rm HI} \rangle = 10^{-4.19}$  $(z=5.4)$',transform=ax[2,1].transAxes, ha='right', va='bottom', fontsize=14)
"""ax[0,0].set_xscale('log')
ax[1,0].set_xscale('log')
ax[2,0].set_xscale('log')
ax[2,1].set_xscale('log')
ax[1,1].set_xscale('log')
ax[0,1].set_xscale('log')"""

#ax[1,1].set_ylim(top=0.0075, bottom=0)
xlim_low = 0.0
xlim_hi = 0.3

ax[0,0].set_xlim(left=xlim_low, right=xlim_hi)
ax[1,0].set_xlim(left=xlim_low, right=xlim_hi)
ax[2,0].set_xlim(left=xlim_low, right=xlim_hi)
ax[2,1].set_xlim(left=xlim_low, right=xlim_hi)
ax[1,1].set_xlim(left=xlim_low, right=xlim_hi)
ax[0,1].set_xlim(left=xlim_low, right=xlim_hi)



plt.tight_layout()
plt.subplots_adjust(left=0.05)

plt.savefig(base + 'fig9.pdf', dpi=300)
plt.show()
"""

fig, ax = plt.subplots(figsize=(fig_width * 2, fig_height * 1.5),

                        nrows=2,

                        ncols=1, sharex=True)
ax[0].plot(mfp_class10_rand.posaxis_p, np.log10(np.mean(mfp_class10_rand.xHI, axis=0)), c='blue', label=r'$\langle x_{HI}\rangle =0.9$')
ax[0].plot(mfp_class6_rand.posaxis_p, np.log10(np.mean(mfp_class6_rand.xHI, axis=0)), c='green', label=r'$\langle x_{HI}\rangle =10^{-1.3}$')
ax[0].plot(mfp_class54_rand.posaxis_p, np.log10(np.mean(mfp_class54_rand.xHI, axis=0)), c='red', label=r'$\langle x_{HI}\rangle =10^{-4.2}$')
ax[0].set_ylabel(r'$\log_{10}\langle  x_{HI} \rangle $',fontsize=16)
ax[1].plot(mfp_class10_rand.posaxis_p, np.mean(np.log10(mfp_class10_rand.xHI), axis=0), c='blue')#, label=r'$\langle x_{HI}\rangle =0.9$')
ax[1].plot(mfp_class6_rand.posaxis_p, np.mean(np.log10(mfp_class6_rand.xHI), axis=0), c='green')#, label=r'$\langle x_{HI}\rangle =10^{-1.3}$')
ax[1].plot(mfp_class54_rand.posaxis_p, np.mean(np.log10(mfp_class54_rand.xHI), axis=0), c='red')
ax[1].set_ylabel(r'$\langle \log_{10} x_{HI} \rangle $',fontsize=16)
ax[1].set_xlabel(r'$R [pMpc]$',fontsize=16)
ax[0].legend(frameon=False, fontsize = 14)
plt.tight_layout()
plt.savefig(base + 'why_doesnt_my_graph_work3.pdf', dpi=300)
plt.show()"""

fig, ax = plt.subplots(figsize=(fig_width*1.1, fig_height * 1.5),

                        nrows=nrows,

                        ncols=1, sharex= 'col')
ax[0].plot(RpMpc_rebin10_3, np.log10(avg_density_rebin10_3), c='cyan', label=r'$\langle M/M_\odot h^{-1} \rangle=10^{9.9}$')
ax[0].plot(RpMpc_rebin10_2, np.log10(avg_density_rebin10_2), c='purple', label=r'$\langle M/M_\odot h^{-1} \rangle=10^{9.4}$')
ax[0].plot(RpMpc_rebin10_1, np.log10(avg_density_rebin10_1), c='red', label=r'$\langle M/M_\odot h^{-1} \rangle=10^{9.1}$')
ax[0].plot(mfp_class10_rand.posaxis_p, np.log10(np.mean(mfp_class10_rand.delta, axis=0)), c='green')
ax[0].set_ylabel(r'$ \log_{10}\langle \Delta \rangle $',fontsize=16) 
ax[0].legend(frameon=False, fontsize = 14)
ax[0].text(0.95,0.2, r'$\langle x_{\rm HI} \rangle = 0.90$',transform=ax[0].transAxes, ha='right', va='bottom', fontsize=14)


ax[1].plot(RpMpc_rebin6_3, np.log10(avg_density_rebin6_3), c='cyan', label=r'$\langle M/M_\odot h^{-1} \rangle=10^{10.7}$')
ax[1].plot(RpMpc_rebin6_2, np.log10(avg_density_rebin6_2), c='purple', label=r'$\langle M/M_\odot h^{-1} \rangle=10^{10.0}$')
ax[1].plot(RpMpc_rebin6_1, np.log10(avg_density_rebin6_1), c='red', label=r'$\langle M/M_\odot h^{-1} \rangle=10^{9.3}$')
ax[1].plot(mfp_class6_rand.posaxis_p, np.log10(np.mean(mfp_class6_rand.delta, axis=0)), c='green')
ax[1].set_ylabel(r'$ \log_{10} \langle\Delta \rangle $',fontsize=16)
ax[1].legend(frameon=False, fontsize = 14)
ax[1].text(0.95,0.2, r'$\langle x_{\rm HI} \rangle = 0.05$',transform=ax[1].transAxes, ha='right', va='bottom', fontsize=14)

ax[2].plot(RpMpc_rebin54_3, np.log10(avg_density_rebin54_3), c='cyan', label=r'$\langle M/M_\odot h^{-1} \rangle=10^{10.9}$')
ax[2].plot(RpMpc_rebin54_2, np.log10(avg_density_rebin54_2), c='purple', label=r'$\langle M/M_\odot h^{-1} \rangle=10^{10.1}$')
ax[2].plot(RpMpc_rebin54_1, np.log10(avg_density_rebin54_1), c='red', label=r'$\langle M/M_\odot h^{-1} \rangle=10^{9.3}$')
ax[2].plot(mfp_class54_rand.posaxis_p, np.log10(np.mean(mfp_class54_rand.delta, axis=0)), c='green')
ax[2].set_ylabel(r'$\log_{10} \langle \Delta \rangle $',fontsize=16)
ax[2].set_xlabel('R [pMpc]',fontsize=16)
ax[2].legend(frameon=False, fontsize = 14)
ax[2].text(0.95,0.2, r'$\langle x_{\rm HI} \rangle = 10^{-4.19}$',transform=ax[2].transAxes, ha='right', va='bottom', fontsize=14)


#ax[1,1].set_ylim(top=0.0075, bottom=0)
xlim_low = 0.0
xlim_hi = 0.3

ax[0].set_xlim(left=xlim_low, right=xlim_hi)
ax[1].set_xlim(left=xlim_low, right=xlim_hi)
ax[2].set_xlim(left=xlim_low, right=xlim_hi)



plt.tight_layout()
#plt.subplots_adjust(left=0.05)

plt.savefig(base + 'density_los_for_talk.pdf', dpi=300)
plt.show()

fig, ax = plt.subplots(figsize=(fig_width, fig_height * 1.5),

                        nrows=nrows,

                        ncols=1, sharex= 'col')

ax[0].plot(RpMpc_rebin10_3, np.log10(avg_H1frac_rebin10_3), c='cyan', label=r'$\langle M/M_\odot h^{-1} \rangle=10^{9.8}$')
ax[0].plot(RpMpc_rebin10_2, np.log10(avg_H1frac_rebin10_2), c='purple', label=r'$\langle M/M_\odot h^{-1} \rangle=10^{9.4}$')
ax[0].plot(RpMpc_rebin10_1, np.log10(avg_H1frac_rebin10_1), c='red', label=r'$\langle M/M_\odot h^{-1} \rangle=10^{9.1}$')
ax[0].plot(mfp_class10_rand.posaxis_p, np.log10(np.mean(mfp_class10_rand.xHI, axis=0)), c='green')
ax[0].set_ylabel(r'$ \log_{10}\langle x_{\rm HI} \rangle $',fontsize=16) 
#ax[0].legend(frameon=False, fontsize = 14)
ax[0].text(0.95,0.45, r'$\langle x_{\rm HI} \rangle = 0.90$',transform=ax[0].transAxes, ha='right', va='bottom', fontsize=14)

ax[1].plot(RpMpc_rebin6_3, np.log10(avg_H1frac_rebin6_3), c='cyan', label=r'$\langle M/M_\odot h^{-1} \rangle=10^{10.7}$')
ax[1].plot(RpMpc_rebin6_2, np.log10(avg_H1frac_rebin6_2), c='purple', label=r'$\langle M/M_\odot h^{-1} \rangle=10^{10.0}$')
ax[1].plot(RpMpc_rebin6_1, np.log10(avg_H1frac_rebin6_1), c='red', label=r'$\langle M/M_\odot h^{-1} \rangle=10^{9.3}$')
ax[1].plot(mfp_class6_rand.posaxis_p, np.log10(np.mean(mfp_class6_rand.xHI, axis=0)), c='green')
ax[1].set_ylabel(r'$\log_{10}\langle  x_{\rm HI} \rangle $',fontsize=16)
#ax[1].legend(frameon=False, fontsize = 14)
ax[1].text(0.95,0.5, r'$\langle x_{\rm HI} \rangle = 0.05$',transform=ax[1].transAxes, ha='right', va='bottom', fontsize=14)

ax[2].plot(RpMpc_rebin54_3, np.log10(avg_H1frac_rebin54_3), c='cyan', label=r'$\langle M/M_\odot h^{-1} \rangle=10^{10.9}$')
ax[2].plot(RpMpc_rebin54_2, np.log10(avg_H1frac_rebin54_2), c='purple', label=r'$\langle M/M_\odot h^{-1} \rangle=10^{10.1}$')
ax[2].plot(RpMpc_rebin54_1, np.log10(avg_H1frac_rebin54_1), c='red', label=r'$\langle M/M_\odot h^{-1} \rangle=10^{9.3}$')
ax[2].plot(mfp_class54_rand.posaxis_p, np.log10(np.mean(mfp_class54_rand.xHI, axis=0)), c='green')
ax[2].set_ylabel(r'$\log_{10}\langle  x_{\rm HI} \rangle $',fontsize=16)
ax[2].set_xlabel('R [pMpc]',fontsize=16)
#ax[2].legend(frameon=False, fontsize = 14)
ax[2].text(0.95,0.45, r'$\langle x_{\rm HI} \rangle = 10^{-4.19}$',transform=ax[2].transAxes, ha='right', va='bottom', fontsize=14)


#ax[1,1].set_ylim(top=0.0075, bottom=0)
xlim_low = 0.0
xlim_hi = 0.3

ax[0].set_xlim(left=xlim_low, right=xlim_hi)
ax[1].set_xlim(left=xlim_low, right=xlim_hi)
ax[2].set_xlim(left=xlim_low, right=xlim_hi)



plt.tight_layout()
#plt.subplots_adjust(left=0.05)

plt.savefig(base + 'xHI_for_talk.pdf', dpi=300)
plt.show()
