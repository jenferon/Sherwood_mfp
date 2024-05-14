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

ncols = 1

#import data for graph

"""high res halo pdf"""
base = '/home/ppxjf3/planck1_40_2048_RTzrfit_los/'
fp_100_z57 = np.loadtxt(base + 'pdf_data/fp_100.txt')
fp_105_z57 = np.loadtxt(base + 'pdf_data/fp_105.txt')
fp_110_z57 = np.loadtxt(base + 'pdf_data/fp_110.txt')
fp_rand_z57 = np.loadtxt(base + 'pdf_data/fp_rand.txt')

pdf_100_z57 ,bins_100_z57, centers_100_z57 = pdf_calc(-6.0, 3, 25, np.log10(fp_100_z57))
pdf_105_z57 ,bins_105_z57, centers_105_z57 = pdf_calc(-6.0, 3, 25, np.log10(fp_105_z57))
pdf_110_z57 ,bins_110_z57, centers_110_z57 = pdf_calc(-6.0, 3, 20, np.log10(fp_110_z57))
pdf_rand_z57, bins_rand_z57, centers_rand_z57 = pdf_calc(-6.0, 3, 20, np.log10(fp_rand_z57))

"""lower res data"""
base = '/home/ppxjf3/planck1_40_512_RTzrfit/pdf_data'
fp_100_z57_lr = np.loadtxt(base + '/fp_100.txt')
fp_105_z57_lr = np.loadtxt(base + '/fp_105.txt')
fp_110_z57_lr = np.loadtxt(base + '/fp_110.txt')
fp_rand_z57_lr = np.loadtxt(base + '/fp_rand.txt')

pdf_100_z57_lr ,bins_100_z57_lr, centers_100_z57_lr = pdf_calc(-6.5, 2, 25, np.log10(fp_100_z57_lr))
pdf_105_z57_lr ,bins_105_z57_lr, centers_105_z57_lr = pdf_calc(-6.5, 2, 25, np.log10(fp_105_z57_lr))
pdf_110_z57_lr ,bins_110_z57_lr, centers_110_z57_lr = pdf_calc(-6.5, 2, 20, np.log10(fp_110_z57_lr))
pdf_rand_z57_lr, bins_rand_z57_lr, centers_rand_z57_lr = pdf_calc(-6.5, 2, 20, np.log10(fp_rand_z57_lr))

fig, ax = plt.subplots(figsize=(fig_width*2, fig_height),nrows=nrows,ncols=ncols)  

#ax.plot(centers_100_z57, pdf_100_z57, color = 'm', label = r'$> 10^{10.0} M_\odot h^{-1}$')
#ax.plot(centers_105_z57, pdf_105_z57, color = 'c', label = r'$> 10^{10.5} M_\odot h^{-1}$')
ax.plot(centers_110_z57, pdf_110_z57, color = 'r', label = r'$> 10^{11.0} M_\odot h^{-1}$')
ax.plot(centers_rand_z57, pdf_rand_z57, color = 'g', label = 'Random LoS')

#ax.scatter(centers_100_z57[np.abs(centers_100_z57 - np.log10(np.mean(fp_100_z57))).argmin()], pdf_100_z57[np.abs(centers_100_z57 - np.log10(np.mean(fp_100_z57))).argmin()], marker = 'X', color = 'm', s = 50)
#ax.scatter(centers_105_z57[np.abs(centers_105_z57 - np.log10(np.mean(fp_105_z57))).argmin()],pdf_105_z57[np.abs(centers_105_z57 - np.log10(np.mean(fp_105_z57))).argmin()], marker = 'X', color = 'c', s = 50)
ax.scatter(centers_110_z57[np.abs(centers_110_z57 - np.log10(np.mean(fp_110_z57))).argmin()],pdf_110_z57[np.abs(centers_110_z57 - np.log10(np.mean(fp_110_z57))).argmin()], marker = 'X', color = 'r', s = 50)
ax.scatter(centers_rand_z57[np.abs(centers_rand_z57 - np.log10(np.mean(fp_rand_z57))).argmin()],pdf_rand_z57[np.abs(centers_rand_z57 - np.log10(np.mean(fp_rand_z57))).argmin()], marker = 'X', color = 'g', s = 50)


#ax.plot(centers_100_z57_lr, pdf_100_z57_lr, color = 'm', linestyle = '--')
#ax.plot(centers_105_z57_lr, pdf_105_z57_lr, color = 'c',  linestyle = '--')
ax.plot(centers_110_z57_lr, pdf_110_z57_lr, color = 'r', linestyle = '--')
ax.plot(centers_rand_z57_lr, pdf_rand_z57_lr, color = 'g', linestyle = '--')

#ax.scatter(centers_100_z57_lr[np.abs(centers_100_z57_lr - np.log10(np.mean(fp_100_z57_lr))).argmin()], pdf_100_z57_lr[np.abs(centers_100_z57_lr - np.log10(np.mean(fp_100_z57_lr))).argmin()], marker = 'o', color = 'm', s = 50)
#ax.scatter(centers_105_z57_lr[np.abs(centers_105_z57_lr - np.log10(np.mean(fp_105_z57_lr))).argmin()],pdf_105_z57_lr[np.abs(centers_105_z57_lr - np.log10(np.mean(fp_105_z57_lr))).argmin()], marker = 'o', color = 'c', s = 50)
ax.scatter(centers_110_z57_lr[np.abs(centers_110_z57_lr - np.log10(np.mean(fp_110_z57_lr))).argmin()],pdf_110_z57_lr[np.abs(centers_110_z57_lr - np.log10(np.mean(fp_110_z57_lr))).argmin()], marker = 'o', color = 'r', s = 50)
ax.scatter(centers_rand_z57_lr[np.abs(centers_rand_z57_lr - np.log10(np.mean(fp_rand_z57_lr))).argmin()],pdf_rand_z57_lr[np.abs(centers_rand_z57_lr - np.log10(np.mean(fp_rand_z57_lr))).argmin()], marker = 'o', color = 'g', s = 50)
ax.set_ylabel(r'$p(log_{10}[\lambda_{fp}$ / pMpc])', fontsize = 16)
ax.set_xlabel(r'$log_{10}[\lambda_{fp}$ / pMpc]', fontsize = 16)
ax.set_yscale('log')
plt.show()