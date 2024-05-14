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
from plotting_dictonary import dictionary

#argparse variables
filename_los = 'filenames.txt'
base1 = 'planck1_40_2048_RTzrfit_los'
base2 = 'planck1_40_2048_RTzr53'
base3 = 'planck1_160_2048_RTzrfit_T40_los'
base4 = 'planck1_40_512_RTzrfit'

#read in text file of file names
lines_los1 = np.genfromtxt(base1 + '/' + filename_los, dtype='str')
lines_los2 = np.genfromtxt(base2 + '/' + filename_los, dtype='str')
lines_los3 = np.genfromtxt(base3 + '/' + filename_los, dtype='str')
lines_los4 = np.genfromtxt(base4 + '/' + filename_los, dtype='str')

#files to store data 
z = np.array([])
av_H1frac1 = np.array([])
av_H1frac2 = np.array([])
av_H1frac3 = np.array([])
av_H1frac4 = np.array([])
for k in range(0, len(lines_los1)):
    # File directory and name for lyman alpha spectra things
    file1  = lines_los1[k]

    # Open the binary file
    readdata = open(base1 + '/' + file1,"rb")

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
    
    if ztime < 13:
        av_H1frac1 = np.append(av_H1frac1, np.mean(H1frac))
        z = np.append(z,ztime)
        # File directory and name for lyman alpha spectra things
        file2  = lines_los2[k]

        # Open the binary file
        readdata = open(base2 + '/' + file1,"rb")

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
        av_H1frac2 = np.append(av_H1frac2, np.mean(H1frac))
        # File directory and name for lyman alpha spectra things
        file3  = lines_los3[k]

        # Open the binary file
        readdata = open(base3 + '/' + file1,"rb")

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
        av_H1frac3 = np.append(av_H1frac3, np.mean(H1frac))
        
        # File directory and name for lyman alpha spectra things
        file4  = lines_los4[k]

        # Open the binary file
        readdata = open(base4 + '/' + file4,"rb")

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
        av_H1frac4 = np.append(av_H1frac4, np.mean(H1frac))
        
    else:
        continue
    

#obvs data 
#Yang et al 2020

x1 = np.array([5.4, 5.6, 5.8, 6.0, 6.2]) 
y1 = np.array([0.00003, 0.000048, 0.000057, 0.000076, 0.000070])
err1 = np.array([(0.0000048,0.000008),(0.0000124, 0.000005),(0.00006518, 0.0000412), (0.00004, 0.000014), (0.0000437, 0.0000079)]).T

#Fan 2006
x2 = np.array([5.03, 5.25, 5.45, 5.65, 5.85, 6.1 ])
y2 = np.array([0.000033, 0.00004, 0.00004, 0.000058, 0.000089, 0.00037])
err2 = np.array([(0.000013,0.0000095), (0.0000178,0.0000153), (0.000021,0.000018), (0.000037, 0.000026), (0.000096, 0.00004), (0.00049, 0.00028)]).T

#McGreer 2015
y_15 = np.array([0.056, 0.368, 0.0382])
x_15 = np.array([5.86, 6.06, 5.58])
yerr_15 = [0.04,0.2,0.05]

#Greig 2018
y_18 = 0.4
x_18 = 7.1
err_18 = np.array([(0.21,0.19)]).T

#Greig 2022
x_22 = np.array([7])
y_22 = np.array([0.64])
err_22 = np.array([(0.19,0.23)]).T

#Jin 2023
x_23 = np.array([5.5, 5.69, 5.89, 6.1, 6.3, 6.49, 6.69])
y_23 = np.array([0.076, 0.16, 0.28, 0.68, 0.78, 0.86, 0.94])
err_23 =  np.array([(0.08,0.08), (0.13,0.14), (0.07, 0.095), (0.055,0.06), (0.05,0.06), (0.03,0.04), (0.06,0.09)]).T

#Zhu 2022
y_z22 = np.array([0.167, 0.281,0.045])
x_z22 = np.array([5.76, 5.94, 5.53])
err_z22 = np.array([0.09, 0.09, 0.03])

#Umeda 2023
x_u23 = np.array([7.14, 7.452, 7.96, 9.801])
y_u23 = np.array([0.46, 0.54, 0.63, 0.83])
yerr_u23 = np.array([(0.32, 0.36),(0.36,0.32),(0.36,0.26),(0.21,0.12)])
xerr_u23 = np.array([(0.076, 0.039), (0.251,0.1),(0.277,0.586),(1.164,1.599)])

np.savetxt('/home/ppxjf3/xh1_data/redshifts.txt', z)
np.savetxt('/home/ppxjf3/xh1_data/zr57_xh1.txt', 1-av_H1frac1)
np.savetxt('/home/ppxjf3/xh1_data/zr57_T40_xh1.txt', 1-av_H1frac2)
np.savetxt('/home/ppxjf3/xh1_data/zr53_xh1.txt', 1-av_H1frac3)
np.savetxt('/home/ppxjf3/xh1_data/zr57_lowres_xh1.txt', 1-av_H1frac4)

"""

fig, ax = plt.subplots(1,1)
obs1 = plt.errorbar(x_18, 1-y_18, yerr= err_18, linestyle ='None', label = 'Greig et al, 2018', capsize=6, color = 'b',  fmt="X")

sim1, = plt.plot(z, 1- av_H1frac1, label=  dictionary['zr57']['label'], linestyle =  dictionary['zr57']['ls'], color =  dictionary['zr57']['c'])
sim2, = plt.plot(z, 1- av_H1frac3, label= dictionary['zr57_T40']['label'], linestyle = dictionary['zr57_T40']['ls'], color = dictionary['zr57_T40']['c'])
sim3, = plt.plot(z, 1- av_H1frac2, label= dictionary['zr53']['label'], linestyle= dictionary['zr53']['ls'] , color = dictionary['zr53']['c'])
sim4, = plt.plot(z, 1- av_H1frac4, label= dictionary['zr57_low_res']['label'], linestyle= dictionary['zr57_low_res']['ls'] , color = dictionary['zr57_low_res']['c'])
#obs2 = plt.errorbar(x1,1-y1,yerr= err1, linestyle ='None', label = 'Yang et al, 2020', capsize=6, color = 'k')
obs3 = plt.errorbar(x_15, 1-y_15, yerr= yerr_15, linestyle ='None', label = 'McGreer+15',   lolims=True, color = dictionary['obs1']['c'], fmt=dictionary['obs1']['fmt'])
#obs4 = plt.errorbar(x2, 1-y2, linestyle ='None', yerr = err2, label = 'Fan et al, 2006', capsize=6, color = 'r')
obs5 = plt.errorbar(x_23, 1- y_23, linestyle ='None', yerr = err_23, label = 'Jin+23', lolims=True, color = dictionary['obs2']['c'], fmt=dictionary['obs2']['fmt'])
obs6 = plt.errorbar(x_22, 1-y_22, yerr= err_22, linestyle ='None', capsize=6, label = 'Greig+22', color = dictionary['obs3']['c'],  fmt=dictionary['obs3']['fmt'])
obs7 = plt.errorbar(x_z22, 1- y_z22, linestyle ='None', yerr = err_z22, label = 'Zhu+22', lolims=True, color = dictionary['obs4']['c'], fmt=dictionary['obs4']['fmt'])
obs8 = plt.errorbar(x_u23, 1- y_u23, linestyle ='None', yerr = yerr_u23, xerr=xerr_u23, label = 'Umeda+23', lolims=True, color = 'mediumspringgreen', fmt='*')
plt.ylabel(r'1 - $\langle x_{H1} \rangle$', fontsize = 16)
plt.xlabel('z', fontsize = 16)
plt.legend(frameon=False, fontsize = 11,  loc='lower left')
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
#plt.yscale('log')
plt.xlim(right=7.5)
plt.savefig('/home/ppxjf3/xh1_data/xh1_comparision.pdf', dpi=300, bbox_inches='tight')
plt.show()


   """