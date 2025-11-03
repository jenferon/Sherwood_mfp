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
parser.add_argument('-t','--filename_tau', help='Files to calculate mfp for', required=True)
parser.add_argument('-l','--filename_los', help='Files to calculate mfp for', required=True)
parser.add_argument('-b','--base', help='base', required=True)
args = parser.parse_args()


#read in text file of file names
lines_tau = np.genfromtxt(str(args.base) + '/' + str(args.filename_tau), dtype='str')
lines_los = np.genfromtxt(str(args.base) + '/' + str(args.filename_los), dtype='str')

#files to store data 
z = np.array([])
F = np.array([])
print(len(lines_los))
print(len(lines_tau))
for k in range(0, len(lines_los)):
    # File directory and name for lyman alpha spectra things
    file1  = lines_los[k]

    # Open the binary file
    readdata = open(str(args.base) + '/' + file1,"rb")

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
    if ztime < 6.3:
        file2  = lines_tau[k]

        # Open the binary file
        readdata = open(str(args.base) + '/' + file2,"rb")
        # H1 Ly-alpha optical depth
        tau_Lya = np.fromfile(readdata,dtype=np.double,count=nbins[0]*numlos[0]) # dimensionless
        readdata.close()
        
        #take 50 cMPc/h steps f
        # or each sightline and then average throughout box 
        step = 50/(box100/1000)
        delta = int(nbins*step)
        print(delta )
        F_temp = np.array([])
        for i in range(0,int((nbins[0]*numlos[0])/delta)):
            F_temp = np.append(F_temp, np.mean(np.exp(-tau_Lya[i*delta : (i+1)*delta])))
            
        F = np.append(F, np.mean(F_temp))
        z = np.append(z,ztime)
    else:
        continue


plt.subplots(1,1)
plt.scatter(z, F)
plt.ylabel(r'$<F>$')
plt.xlabel('z')
#plt.legend()
#plt.savefig('/home/ppxjf3/mfp_data/' + str(args.savefile)+'.pdf')
plt.show()
"""
np.savetxt('/home/ppxjf3/optical_depths/' + str(args.base) +'.txt', (z,F))
for i in range(0,len(F)):
    print('at z = ' + str(z[i]) + ' <F> = ' +str(F[i]))

  """     