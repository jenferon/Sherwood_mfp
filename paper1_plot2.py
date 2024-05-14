#make first paper overview plot

import numpy as np
from plotting_dictonary import dictionary
import matplotlib.pyplot as plt
import matplotlib

#Paper specific Matploblib settings
import matplotlib.font_manager
import matplotlib as mpl




def pdf_calc (bin_min, bin_max, nbins, data):
	step = (bin_max - bin_min)/nbins
	#print(step)
	bins = np.arange(bin_min, bin_max, step)
	N = np.array([])
	for j in range (0, nbins-1):
		
		Ni =np.where((data>= bins[j]) & (data < bins[j+1]))
		N = np.append(N,len(Ni[0]))
	
	N_tot = np.sum(N)
	pdf = N/(N_tot*step)
	sum = 0
	for k in range(0,len(pdf)):
		sum += pdf[k]*(bins[k+1]-bins[k])
	#print(sum)
	centers=0.5*(bins[:-1]+bins[1:])
	return pdf, bins, centers

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

#DPI of MNRAS is 300
mpl.rcParams['figure.dpi'] = 300/scale
# Set up a figure with three panels, with one rows and three columns

nrows = 1

ncols = 2


#import data for graph

"""halo pdf"""
base = '/home/ppxjf3/planck1_40_2048_RTzrfit_los/'
fp_100 = np.loadtxt(base + 'pdf_data/fp_100.txt')
fp_105 = np.loadtxt(base + 'pdf_data/fp_105.txt')
fp_110 = np.loadtxt(base + 'pdf_data/fp_110.txt')
fp_rand = np.loadtxt(base + 'pdf_data/fp_rand.txt')

pdf_100 ,bins_100, centers_100 = pdf_calc(-6.0, 3, 25, np.log10(fp_100))
pdf_105 ,bins_105, centers_105 = pdf_calc(-6.0, 3, 25, np.log10(fp_105))
pdf_110 ,bins_110, centers_110 = pdf_calc(-6.0, 3, 25, np.log10(fp_110))
pdf_rand, bins_rand, centers_rand = pdf_calc(-6.0, 3, 25, np.log10(fp_rand))

"""Redshift evolution graph"""
fp_z10 = np.loadtxt(base+'pdf_data/fp_z=10.0_rand.txt')
fp_z8 = np.loadtxt(base+'pdf_data/fp_z=8.0_rand.txt')
fp_z7 = np.loadtxt(base+'pdf_data/fp_z=7.0_rand.txt')
fp_z6 = np.loadtxt(base+'pdf_data/fp_z=6.0_rand.txt')
fp_z54 = np.loadtxt(base+'pdf_data/fp_z=5.4_rand.txt')

pdf_z10 ,bins_z10, centers_z10 = pdf_calc(-6.0, 3, 25, np.log10(fp_z10))
pdf_z8 ,bins_z8, centers_z8 = pdf_calc(-6.0, 3, 25, np.log10(fp_z8))
pdf_z7 ,bins_z7, centers_z7 = pdf_calc(-6.0, 3, 25, np.log10(fp_z7))
pdf_z6, bins_z6, centers_z6 = pdf_calc(-6.0, 3, 25, np.log10(fp_z6))
pdf_z54, bins_z54, centers_z54 = pdf_calc(-6.0, 3, 25, np.log10(fp_z54))


"""Los averages"""
delta_100 = np.loadtxt(base+'/los_data/delta_100.txt')
delta_105 = np.loadtxt(base+'/los_data/delta_105.txt')
delta_110 = np.loadtxt(base+'/los_data/delta_110.txt')
delta_rand = np.loadtxt(base+'/los_data/delta_rand.txt')
xHI_100 = np.loadtxt(base+'/los_data/xHI_100.txt')
xHI_105 = np.loadtxt(base+'/los_data/xHI_105.txt')
xHI_110 = np.loadtxt(base+'/los_data/xHI_110.txt')
xHI_rand = np.loadtxt(base+'/los_data/xHI_rand.txt')
posaxis_100 = np.loadtxt(base+'/los_data/posaxis_100.txt')
posaxis_105 = np.loadtxt(base+'/los_data/posaxis_105.txt')
posaxis_110 = np.loadtxt(base+'/los_data/posaxis_110.txt')
posaxis_rand = np.loadtxt(base+'/los_data/posaxis_rand.txt')

# axs is a numpy array with dimension (nrows, ncols)

fig, ax = plt.subplots(figsize=(fig_width*2, fig_height),

                        nrows=nrows,

                        ncols=ncols)    

ax[0].plot(centers_100, pdf_100, color = 'm', label = r'$> 10^{10.0} M_\odot h^{-1}$')
ax[0].plot(centers_105, pdf_105, color = 'c', label = r'$> 10^{10.5} M_\odot h^{-1}$')
ax[0].plot(centers_110, pdf_110, color = 'r', label = r'$> 10^{11.0} M_\odot h^{-1}$')
ax[0].plot(centers_rand, pdf_rand, color = 'g', label = 'Random LoS')

ax[0].scatter(centers_100[np.abs(centers_100 - np.log10(np.mean(fp_100))).argmin()], pdf_100[np.abs(centers_100 - np.log10(np.mean(fp_100))).argmin()], marker = 'X', color = 'm', s = 50)
ax[0].scatter(centers_105[np.abs(centers_105 - np.log10(np.mean(fp_105))).argmin()],pdf_105[np.abs(centers_105 - np.log10(np.mean(fp_105))).argmin()], marker = 'X', color = 'c', s = 50)
ax[0].scatter(centers_110[np.abs(centers_110 - np.log10(np.mean(fp_110))).argmin()],pdf_110[np.abs(centers_110 - np.log10(np.mean(fp_110))).argmin()], marker = 'X', color = 'r', s = 50)
ax[0].scatter(centers_rand[np.abs(centers_rand - np.log10(np.mean(fp_rand))).argmin()],pdf_rand[np.abs(centers_rand - np.log10(np.mean(fp_rand))).argmin()], marker = 'X', color = 'g', s = 50)

ax[0].set_ylabel(r'$p(log_{10}[\lambda_{fp}$ / pMpc])', fontsize = 16)
ax[0].set_xlabel(r'$log_{10}[\lambda_{fp}$ / pMpc]', fontsize = 16)
ax[0].set_yscale('log')


"""colour = ['b', 'm', 'darkorange', 'g', 'r']
ax[1].plot(centers_z10, pdf_z10, color = colour[0], label = r'z=10')
ax[1].plot(centers_z8, pdf_z8, color = colour[1], label = r'z=8')
ax[1].plot(centers_z7, pdf_z7, color = colour[2], label = r'z=7')
ax[1].plot(centers_z6, pdf_z6, color = colour[3], label = r'z=6')
ax[1].plot(centers_z54, pdf_z54, color = colour[4], label = r'z=5.4')

ax[1].scatter(centers_z10[np.abs(centers_z10 - np.log10(np.mean(fp_z10))).argmin()], pdf_z10[np.abs(centers_z10 - np.log10(np.mean(fp_z10))).argmin()], marker = 'X', color = colour[0], s = 50)
ax[1].scatter(centers_z8[np.abs(centers_z8 - np.log10(np.mean(fp_z8))).argmin()],pdf_z8[np.abs(centers_z8 - np.log10(np.mean(fp_z8))).argmin()], marker = 'X', color =  colour[1], s = 50)
ax[1].scatter(centers_z7[np.abs(centers_z7 - np.log10(np.mean(fp_z7))).argmin()],pdf_z7[np.abs(centers_z7 - np.log10(np.mean(fp_z7))).argmin()], marker = 'X', color =  colour[2], s = 50)
ax[1].scatter(centers_z6[np.abs(centers_z6 - np.log10(np.mean(fp_z6))).argmin()],pdf_z6[np.abs(centers_z6 - np.log10(np.mean(fp_z6))).argmin()], marker = 'X', color =  colour[3], s = 50)
ax[1].scatter(centers_z54[np.abs(centers_z54 - np.log10(np.mean(fp_z54))).argmin()], pdf_z54[np.abs(centers_z54 - np.log10(np.mean(fp_z54))).argmin()], marker = 'X', color =  colour[4], s = 50)
ax[1].set_ylabel(r'$p(log_{10}[\lambda_{fp}$ / pMpc])', fontsize = 16)
ax[1].set_xlabel(r'$log_{10}[\lambda_{fp}$ / pMpc]', fontsize = 16)
ax[1].set_yscale('log')
ax[1].legend(frameon=False, fontsize = 11)"""

ax[1].plot(posaxis_rand, delta_rand ,c='g', label = 'random')
ax[1].plot(posaxis_100, delta_100, c= 'm', label = r'$> 10^{10.0} M_\odot h^{-1}$')
ax[1].plot(posaxis_110, delta_110, c='r', label = r'$> 10^{11.0} M_\odot h^{-1}$')
ax[1].plot(posaxis_105, delta_105, c='c', label = r'$> 10^{10.5} M_\odot h^{-1}$')
ax[1].set_ylabel(r'$\Delta$',fontsize=12) 
ax[1].legend(frameon=False, fontsize = 11)
#ax[0,2].plot(posaxis_rand, xHI_rand , c='g', label = 'random', linestyle = '--')
#ax[0,2].plot(posaxis_100, xHI_100, c= 'm', label = r'$> 10^{10.0} M_\odot h^{-1}$', linestyle = '--')
#ax[0,2].plot(posaxis_110, xHI_110, c='r', label = r'$> 10^{11.0} M_\odot h^{-1}$', linestyle = '--')
#ax[0,2].plot(posaxis_105, xHI_105, c='c', label = r'$> 10^{10.5} M_\odot h^{-1}$', linestyle = '--')
#ax[0,2].set_ylabel(r'$x_{H1} $',fontsize=12)
ax[1].set_yscale('log')
#ax[0,2].legend(frameon=False)
ax[1].set_xlabel(r'$\rm x\,[ckPc\, h^{-1}]$',fontsize=12)
ax[1].set_xlim(left = 0 , right=10000)#

plt.xticks(fontsize=16)
plt.yticks(fontsize=16)

plt.tight_layout()
plt.savefig(base + 'paper1_fig2.pdf', dpi=300, bbox_inches='tight')
#plt.show()