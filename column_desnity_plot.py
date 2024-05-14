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

#load in data for threshold graph
base = '/home/ppxjf3/planck1_40_2048_RTzrfit_los/column_desnity_data/'
NHI_100 = np.loadtxt(base+'NHI_100.txt')
dr_100 = np.loadtxt(base+'dr_100.txt')
NHI_80 = np.loadtxt(base+'NHI_80.txt')
dr_80= np.loadtxt(base+'dr_80.txt')
NHI_70 = np.loadtxt(base+'NHI_70.txt')
dr_70 = np.loadtxt(base+'dr_70.txt')
NHI_60 = np.loadtxt(base+'NHI_60.txt')
dr_60 = np.loadtxt(base+'dr_60.txt')
NHI_54 = np.loadtxt(base+'NHI_54.txt')
dr_54 = np.loadtxt(base+'dr_54.txt')

NHI_100_thresh = np.loadtxt(base+'NHI_100_thresh_0.5.txt')
NHI_80_thresh = np.loadtxt(base+'NHI_80_thresh_0.5.txt')
NHI_70_thresh = np.loadtxt(base+'NHI_70_thresh_0.5.txt')
NHI_60_thresh = np.loadtxt(base+'NHI_60_thresh_0.5.txt')
NHI_54_thresh = np.loadtxt(base+'NHI_54_thresh_0.5.txt')

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

#for window size variation plot
NHI_100_50 = np.loadtxt(base+'NHI_100_w=50.txt')
NHI_100_25 = np.loadtxt(base+'NHI_100_w=25.txt')
NHI_100_100 = np.loadtxt(base+'NHI_100_w=100.txt')
NHI_60_50 = np.loadtxt(base+'NHI_60_w=50.txt')
NHI_60_25 = np.loadtxt(base+'NHI_60_w=25.txt')
NHI_60_100 = np.loadtxt(base+'NHI_60_w=100.txt')
dr_100_50 = np.loadtxt(base+'dr_100_w=50.txt')
dr_100_25 = np.loadtxt(base+'dr_100_w=25.txt')
dr_100_100 = np.loadtxt(base+'dr_100_w=100.txt')
dr_60_50 = np.loadtxt(base+'dr_60_w=50.txt')
dr_60_25 = np.loadtxt(base+'dr_60_w=25.txt')
dr_60_100 = np.loadtxt(base+'dr_60_w=100.txt')

pdf_100_50, bins_100_50, centers_100_50 = pdf_calc(np.min(NHI_100_50), np.max(NHI_100_50), 18, NHI_100_50, dr_100_50)
pdf_100_25, bins_100_25, centers_100_25 = pdf_calc(np.min(NHI_100_25), np.max(NHI_100_25), 18, NHI_100_25, dr_100_25)
pdf_100_100, bins_100_100, centers_100_100 = pdf_calc(np.min(NHI_100_100), np.max(NHI_100_100), 18, NHI_100_100, dr_100_100)
pdf_60_50, bins_60_50, centers_60_50 = pdf_calc(np.min(NHI_60_50), np.max(NHI_60_50), 18, NHI_60_50, dr_60_50)
pdf_60_25, bins_60_25, centers_60_25 = pdf_calc(np.min(NHI_60_25), np.max(NHI_60_25), 18, NHI_60_25, dr_60_25)
pdf_60_100, bins_60_100, centers_60_100 = pdf_calc(np.min(NHI_60_100), np.max(NHI_60_100), 18, NHI_60_100, dr_60_100)

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
fig, ax = plt.subplots(figsize=(fig_width, fig_height),nrows=nrows,ncols=ncols)


ax.plot(centers_100, np.log10(pdf_100), color = colour[0], label=r'$\langle x_{\mathrm{H1}} \rangle = 0.9$')
ax.plot(centers_80, np.log10(pdf_80), color = colour[1], label=r'$\langle x_{\mathrm{H1}} \rangle = 10^{-0.2}$')
ax.plot(centers_70, np.log10(pdf_70), color = colour[2], label=r'$\langle x_{\mathrm{H1}} \rangle = 10^{-0.4}$')
ax.plot(centers_60, np.log10(pdf_60), color = colour[3], label=r'$\langle x_{\mathrm{H1}} \rangle = 10^{-1.3}$')
ax.plot(centers_54, np.log10(pdf_54), color = colour[4], label=r'$\langle x_{\mathrm{H1}} \rangle = 10^{-4.2}$')

ax.plot(centers_54_thresh, np.log10(pdf_54_thresh), linestyle = '--', color = colour[4])#, label=r'$\langle x_{\mathrm{H1}} \rangle = 10^{-4.2}$')#.format(np.log10(av_xHI_54)))
ax.plot(centers_60_thresh, np.log10(pdf_60_thresh), linestyle='--', color = colour[3])#, label=r'$\langle x_{\mathrm{H1}} \rangle = 10^{-1.3}$')#label='xHI=1e{:.2f}'.format(np.log10(av_xHI_60)))
ax.plot(centers_70_thresh, np.log10(pdf_70_thresh), color = colour[2], linestyle='--')#, label=r'$\langle x_{\mathrm{H1}} \rangle = 10^{-0.4}$')#label='xHI=1e{:.2f}'.format(np.log10(av_xHI_70)))
ax.plot(centers_80_thresh, np.log10(pdf_80_thresh),linestyle = '--', color = colour[1])#, label=r'$\langle x_{\mathrm{H1}} \rangle = 10^{-0.2}$')#label='xHI=1e{:.2f}'.format(np.log10(av_xHI_80)))
ax.plot(centers_100_thresh, np.log10(pdf_100_thresh), color = colour[0], linestyle = '--')# label=r'$\langle x_{\mathrm{H1}} \rangle = 0.9$')

ax.set_xlabel(r'$log_{10}[N_{\mathrm{H1}} / \mathrm{cm}^{-2}]$', fontsize = 18)
ax.set_ylabel(r'$log_{10}[\frac{ \partial ^2n}{ \partial log_{10}[N_{\mathrm{H1}}] \partial R} /\mathrm{ pMPc}^{-1}]$', fontsize = 18)
#ax[0].set_xticks(fontsize=16)
#ax[0].set_yticks(fontsize=16)
ax.legend(frameon=False, fontsize = 11)

"""ax[1].plot(centers_100_50, np.log10(pdf_100_50), color = colour[0], label=r'$\langle x_{\mathrm{H1}} \rangle = 0.9$')#label='xHI=1e{:.2f}'.format(np.log10(av_xHI_100)))
ax[1].plot(centers_100_100, np.log10(pdf_100_100), color = colour[0], linestyle = '--')
ax[1].plot(centers_100_25, np.log10(pdf_100_25), color = colour[0], linestyle = ':')
ax[1].plot(centers_60_50, np.log10(pdf_60_50), color = colour[3], label=r'$\langle x_{\mathrm{H1}} \rangle = 10^{-1.3}$')
ax[1].plot(centers_60_100, np.log10(pdf_60_100), color = colour[3], linestyle = '--')
ax[1].plot(centers_60_25, np.log10(pdf_60_25), color = colour[3], linestyle=':')

ax[1].set_xlabel(r'$log_{10}[N_{\mathrm{H1}} /\mathrm{cm}^{-2}]$', fontsize = 18)
ax[1].set_ylabel(r'$log_{10}[\frac{ \partial ^2n}{ \partial log_{10}[N_{\mathrm{H1}}] \partial R} / \mathrm{pMPc}^{-1}]$', fontsize = 18)
#ax[1].set_xticks(fontsize=16)
#ax[1].set_yticks(fontsize=16)
ax[1].legend(frameon=False, fontsize = 11)"""
plt.tight_layout()
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)


plt.savefig(base + 'paper1_fig3.pdf', bbox_inches='tight', pad_inches=0.02)  # remove whitespace
#plt.show()