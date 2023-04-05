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


# File directory and name for lyman alpha spectra things
base = '/home/ppxjf3/planck1_40_2048_RTzrfit_los'
file1  = 'los2048_n5000_z6.400.dat'
# Open the binary file
readdata = open(str(base) + '/' + file1,"rb")

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
readdata.close()
        
file2  = 'tau2048_n5000_z6.400.dat'

# Open the binary file
readdata = open(str(base) + '/' + file2,"rb")
# H1 Ly-alpha optical depth
tau_Lya = np.fromfile(readdata,dtype=np.double,count=nbins[0]*numlos[0]) # dimensionless
print(tau_Lya)
readdata.close()
""" """    
#take 50 cMPc/h steps f
# or each sightline and then average throughout box 
step = 50/(box100/1000)
delta = int(nbins*step)
t_eff = np.array([])
temp = np.array([])
iterations = 1000
for j in range (0, iterations):
            
    index = random.randint(0,len(tau_Lya)-1)
    print(index)
    for num in range(0,delta):
        #print('num = ' + str(num) + ' index = ' + str(index))
        if index == len(tau_Lya):
            #incase we reach the end of the box
            print('end of box')
            index  =0
        temp = np.append(temp, tau_Lya[index])
                

        index += 1     
    t_eff = np.append(t_eff, np.mean(temp))
    print(t_eff[j])
    print(j)

x = np.sort(t_eff)
#print(t_eff.shape)
# get the cdf values of y
y = np.arange(iterations) / iterations
# plotting
plt.xlabel(r'$<\tau_eff>')
plt.ylabel('CDF')
  
plt.title('CDF for 40_2048_RTzrfit')
  
plt.plot(x, y, label = 'z = ' + str(ztime))
plt.legend()
plt.savefig('cdf.pdf')
plt.show()