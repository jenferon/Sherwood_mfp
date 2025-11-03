import numpy as np
import matplotlib.pyplot as plt
import math
from scipy import interpolate
#import seaborn as sns
from mpl_toolkits.axes_grid1 import make_axes_locatable
from astropy import constants as const
import random
import scipy
import scipy.stats
import argparse


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

parser = argparse.ArgumentParser(description='Files and parameters to calculate mfp.')
parser.add_argument('-f','--folder', help='folder with relevant files', required=True)
parser.add_argument('-z','--redshifts', nargs='+', help='list of redshifts to plot', required=True)
args = parser.parse_args()



base = args.folder
z = args.redshifts
itter = 10000
mfp_sher = np.empty([len(z),itter])
xH1 = np.empty([len(z)])

p = 0
for i in z:

    file = base + '/los2048_n5000_z' + str(i) + '00.dat'

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
    -----------------------------------------------------------------------------------------------------------------------------------------------
    Mean free path calculation
    -----------------------------------------------------------------------------------------------------------------------------------------------
    """

        
    #define constants 
    n = number_density_H(h100,density, omegab, Xh, ztime)*1e-6 #cm-3
    n_HI = n*H1frac #cm-3
            
    #define cross section
    cross_sec = 6.3*10**-18 #cm2 for 912 amstrongs

    #empty array for mfp
    mfp= np.array([])

    #convert units 
    cm_pos = (posaxis*kPC*1e2)/((1+ztime)*h100)#proper cm
    delta_r = cm_pos[3] - cm_pos[2]
    
	#random number seeded so the same random numbers are generated each time
    np.random.seed(0) 
    start_value = np.random.randint(len(n_HI)-1, size = itter)


    for j in start_value:
        ii = j
        #print(i)
        tau = 0
        n= 0
        while tau <1:
            tau += cross_sec*n_HI[ii]*delta_r #dimensionless
            #print(tau)
            n += 1
            ii += 1
            if ii == len(n_HI):
            #incase we reach the end of the boxx
                print('end of box')
                ii =0
        mfp =np.append(mfp,(n*delta_r)/tau)
        #print(n)
        #print(j)

    #convert to pMpc from comoving cm 
    proper_mfp = mfp*1e-2 /(MPC)
    mfp_sher[p,:] = proper_mfp
    av_xhI = np.mean(H1frac)
    xH1[p] = np.round(np.log10(av_xhI), decimals=1)
    p +=1
    np.savetxt(base+'/pdf_data/fp_z=' + str(i) +'_rand.txt', proper_mfp)

log_mfp = np.log10(mfp_sher)
#plot pdf of mfp distribution
bin_min = -6.0
bin_max = 3.0
nbins = 25
step = (bin_max - bin_min)/nbins
print(step)
bins = np.arange(bin_min, bin_max, step)
colour = ['b', 'm', 'darkorange', 'g', 'r']
plt.subplots(1,1)
#print(len(z))
for i in range(0,len(z)):    
    #hist, bin_edges  = np.histogram(mfp_sher[i,:], bins = bins)
    N = np.array([])
    for j in range (0,nbins-1):
        
        Ni =np.where((log_mfp[i,:]>= bins[j]) & (log_mfp[i,:] < bins[j+1]))
        #print(np.log10(mfp_sher[i,Ni[0]]))
        N = np.append(N,len(Ni[0]))
    #print(np.log10(mfp_sher[i,:]))
    centers=0.5*(bins[:-1]+bins[1:])
    N_tot = np.sum(N)
    pdf = N/(N_tot*step)
    sum = 0
    for k in range(0,len(pdf)):
        sum += pdf[k]*(bins[k+1]-bins[k])
    print(sum)
    #plt.plot(centers, pdf, linewidth=2,  label='z= '+z[i])
    #hist, bin_edges  = np.histogram(np.log10(mfp_sher[i,:]), bins = bins, density=True)
    #sum = 0
    #for k in range(0,len(hist)):
        #sum += hist[k]*(bin_edges[k+1]-bin_edges[k])
    #print(sum)
    if i == 0:
        plt.plot(centers, pdf, label=r'$\langle x_{\mathrm{H1}} \rangle = 0.9$', c = colour[i])
    if i == 1:
        plt.plot(centers, pdf, label=r'$\langle x_{\mathrm{H1}} \rangle = 10^{-0.2}$', c = colour[i])
    if i == 2:
        plt.plot(centers, pdf, label=r'$\langle x_{\mathrm{H1}} \rangle = 10^{-0.4}$', c = colour[i])
    if i == 3:
        plt.plot(centers, pdf, label=r'$\langle x_{\mathrm{H1}} \rangle = 10^{-1.3}$', c = colour[i])
    if i == 4:
        plt.plot(centers, pdf, label=r'$\langle x_{\mathrm{H1}} \rangle = 10^{-4.2}$', c = colour[i])
    mfp = np.log10(np.mean(mfp_sher[i,:]))
    print('mfp =' +str(mfp))
    
    plt.scatter(centers[np.abs(centers - np.log10(np.mean(mfp_sher[i,:]))).argmin()],pdf[np.abs(centers - np.log10(np.mean(mfp_sher[i,:]))).argmin()], marker = 'X', c = colour[i], s = 50)
plt.ylabel(r'$p(log_{10}[\lambda_{\mathrm{fp}}$ / pMpc])', fontsize = 16)
plt.xlabel(r'$log_{10}[\lambda_{\mathrm{fp}}$ / pMpc]', fontsize = 16)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.yscale('log')
plt.ylim(bottom = 0.001)
plt.xlim(left = -6, right = 2)
plt.legend(frameon=False, fontsize = 8, loc = 'upper left')
plt.savefig('mfp_dist_'+str(base)+'.pdf', dpi=300, bbox_inches='tight')
plt.show()

#np.savetxt('40_512_RTzrfit_6.0_logmfp.txt', log_mfp)

