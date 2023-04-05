import numpy as np
import matplotlib.pyplot as plt
import math
from scipy import interpolate
#import seaborn as sns
from mpl_toolkits.axes_grid1 import make_axes_locatable
from astropy import constants as const
import random
import scipy
import argparse

#argparse variables
parser = argparse.ArgumentParser(description='Files and parameters to calculate mfp.')
parser.add_argument('-f','--folder', help='folder with relevant files', required=True)
parser.add_argument('-N','--itterations', help='Number of itterations in mfp calculation', required=True)
parser.add_argument('-z_up','--high_z', help='upper z limit', required=True)
args = parser.parse_args()

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

    

#read in text file of file names
lines = np.genfromtxt(str(args.folder)+'/filenames.txt', dtype='str')#'filenames.txt',dtype='str')

z = np.array([])
mfp_sher = np.array([])
error = np.array([])


for k in range(0, len(lines)):
    """_
    ------------------------------------------------------------------------------------------------------------------------
    Read in binary file data 
    ------------------------------------------------------------------------------------------------------------------------
    """
    # File directory and name for lyman alpha spectra things
    file  = lines[k]

    # Open the binary file
    readdata = open(str(args.folder) + '/' +file,"rb")
    
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
    print(ztime)
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

    
    if ztime < float(args.high_z):
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


        for j in range (0, int(args.itterations)):
            i = random.randint(0,len(n_HI)-1)
            tau = 0
            n= 0
            while tau <1:
                tau += cross_sec*n_HI[i]*delta_r #dimensionless
                #print(tau)
                n += 1
                i += 1
                if i == len(n_HI):
                    #incase we reach the end of the box
                    i =0
            mfp =np.append(mfp,n*delta_r)
            print(j)

        #convert to pMpc from comoving cm 
        proper_mfp = mfp*1e-2 /(MPC)
        
        #calculate error with the bootstrapping method
        #err = scipy.stats.bootstrap([proper_mfp], np.std, confidence_level=0.9).standard_error
            
        mfp_sher = np.append(mfp_sher, np.mean(proper_mfp))
        z = np.append(z,ztime)
        #error = np.append(error, err)



np.savetxt('/home/ppxjf3/mfp_data/' + str(args.folder)+ '_mfp_N' +str(args.itterations) +'.txt', (z,mfp_sher))
#define points from literature - in pMpc
obs_mfp = np.array([[6, 5.1, 5.16, 4.86, 4.56],[0.75, 9.09, 10.3, 15.1, 22.2], [0.45, 1.28, 1.6, 1.8, 2.3], [0.65, 1.62, 1.6, 1.8, 2.3]])

plt.subplots(1,1)
plt.errorbar(obs_mfp[0,:], obs_mfp[1,:],yerr=obs_mfp[2:4,:], fmt=".", label = 'Observations')
plt.errorbar(z, mfp_sher,  yerr=0 , fmt=".", label = 'Sherwood', color='k')
plt.yscale('log')
plt.ylabel(r'$\lambda_{mfp}$ (pMpc)')
plt.xlabel('z')
plt.legend()
plt.savefig('/home/ppxjf3/mfp_data/'  + str(args.folder)+ '_mfp_N' +str(args.itterations) +'.pdf')
plt.show()

