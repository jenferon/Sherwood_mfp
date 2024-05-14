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
av_gamma1 = np.array([])
av_gamma2 = np.array([])
av_gamma3 = np.array([])
av_gamma4 = np.array([])
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

    Gammaker = np.fromfile(readdata,dtype=np.double,count=nbins[0]*numlos[0]) #/* Gamma_HIc [s^-1] */
    # Close the binary file
    readdata.close()
    
    if ztime < 8:
        av_gamma1 = np.append(av_gamma1, np.mean(Gammaker)*(10**12))
        z = np.append(z,ztime)
        # File directory and name for lyman alpha spectra things
        file2  = lines_los2[k]

        # Open the binary file
        readdata = open(base2 + '/' + file2,"rb")

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
        
        Gammaker = np.fromfile(readdata,dtype=np.double,count=nbins[0]*numlos[0]) #/* Gamma_HIc [s^-1] */

        # Close the binary file
        readdata.close()
        
        av_gamma2 = np.append(av_gamma2, np.mean(Gammaker)*(10**12))
        
        file3  = lines_los3[k]

        # Open the binary file
        readdata = open(base3 + '/' + file3,"rb")

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
        
        Gammaker = np.fromfile(readdata,dtype=np.double,count=nbins[0]*numlos[0]) #/* Gamma_HIc [s^-1] */

        # Close the binary file
        readdata.close()
        av_gamma3 = np.append(av_gamma3, np.mean(Gammaker)*(10**12))
        
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
        
        Gammaker = np.fromfile(readdata,dtype=np.double,count=nbins[0]*numlos[0]) #/* Gamma_HIc [s^-1] */

        # Close the binary file
        readdata.close()
        av_gamma4 = np.append(av_gamma4, np.mean(Gammaker)*(10**12))
        
    else:
        continue
    

#obvs data 
#Gaikward 2023
x_g23 = [4.9, 5.0, 5.1, 5.2, 5.3, 5.4, 5.5, 5.6, 5.7, 5.8, 5.9, 6.0]
y_g23 = [0.501, 0.557, 0.508, 0.502, 0.404, 0.372, 0.344, 0.319, 0.224, 0.178, 0.151, 0.145]
xerr_g23 = 0.05
yerr_g23 = np.array([(0.232, 0.275), (0.218, 0.376), (0.192, 0.324), (0.193, 0.292), (0.147, 0.272), (0.126, 0.217), (0.13, 0.219), (0.12, 0.194), (0.112, 0.223), (0.078, 0.194),(0.079, 0.151), (0.087, 0.157)]).T

#Calverley 2011
y_c11 = [10**-0.15, 10**-0.77] #[0.08, 0.08, 0.03, 0.10, 0.012, 0.23, 0.43, 0.14, 0.2, 1.82, 0.48, 0.95, 1.02, 3.89, 1.88]
x_c11 = [5.04, 6.09] #[6.4189, 6.308, 6.247, 6.02, 6.016, 5.82, 5.81, 5.41, 5.33, 5.2, 5.09, 4.967, 5.886, 4.876, 4.588]
#yerr_c11 = np.array([(0.06, 0.19), (0.05, 0.15), (0.02, 0.1), (0.06, 0.17), (0.08, 0.30), (0.15,0.4), (0.26, 0.65), (0.09, 0.26), (0.13, 0.4), (1.69, 23.6), (0.38, 1.88), (0.52, 1.13), (0.72, 2.51), (2.86, 10.8), (1.33, 4.58)]).T
yerr_c11 = np.array([(0.49,1.023), (0.112,0.257)])
xerr_c11 = np.array([(0.44,0.36), (0.29,0.31)])

#Becker 2013
x_b13 = [4.0, 4.4, 4.75]
y_b13 = [-0.072, -0.019, -0.029]
yerr_b13 = np.array([(0.117,0.135), (0.112,0.14), (0.147,0.156)])

np.savetxt('/home/ppxjf3/gamma_data/gamma_redshifts.txt', z)
np.savetxt('/home/ppxjf3/gamma_data/zr57_xh1.txt', av_gamma1)
np.savetxt('/home/ppxjf3/gamma_data/zr57_T40_xh1.txt', av_gamma2)
np.savetxt('/home/ppxjf3/gamma_data/zr53_xh1.txt', av_gamma3)
np.savetxt('/home/ppxjf3/gamma_data/zr57_lowres_xh1.txt', av_gamma4)
"""
fig, ax = plt.subplots(1,1)


sim1, = plt.plot(z, av_gamma1, label=  dictionary['zr57']['label'], linestyle =  dictionary['zr57']['ls'], color =  dictionary['zr57']['c'])
sim2, = plt.plot(z, av_gamma3, label= dictionary['zr57_T40']['label'], linestyle = dictionary['zr57_T40']['ls'], color = dictionary['zr57_T40']['c'])
sim3, = plt.plot(z, av_gamma2, label= dictionary['zr53']['label'], linestyle= dictionary['zr53']['ls'] , color = dictionary['zr53']['c'])
sim4, = plt.plot(z, av_gamma4, label= dictionary['zr57_low_res']['label'], linestyle= dictionary['zr57_low_res']['ls'] , color = dictionary['zr57_low_res']['c'])
obs1 = plt.errorbar(x_g23, y_g23, yerr= yerr_g23, xerr=xerr_g23, linestyle ='None', label = 'Gaikwad+23', capsize=6, color = dictionary['obs1']['c'],  fmt= dictionary['obs1']['fmt'])
obs2 = plt.errorbar(x_b13, y_b13, yerr= yerr_b13,  linestyle ='None', label = 'Becker+13', capsize=6, color = dictionary['obs2']['c'],  fmt= dictionary['obs2']['fmt'])
obs3 = plt.errorbar(x_c11, y_c11, yerr= yerr_c11, xerr = xerr_c11, linestyle ='None', label = 'Calverley+11', capsize=6, color = dictionary['obs4']['c'],  fmt= dictionary['obs4']['fmt'])

plt.ylabel(r'$\langle \Gamma_{HI} \rangle \times 10^{-12} ($ s$^{-1})$' , fontsize = 16)
plt.xlabel('z', fontsize = 16)
first_legend = ax.legend(handles=[obs1, obs2, obs3], loc='upper right', frameon=False, fontsize = 11)
# Add the legend manually to the Axes.
ax.add_artist(first_legend)
#ax.legend(handles=[sim_1, sim_3, sim_4], loc='lower left', frameon=False, fontsize = 11)
ax.legend(handles=[sim1, sim2, sim3, sim4], loc='lower left', frameon=False, fontsize = 11)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
#plt.yscale('log')
plt.xlim(right=7.5)
plt.yscale('log')
plt.savefig('/home/ppxjf3/mfp_data/gamma_comparision.pdf', dpi=300, bbox_inches='tight')
plt.show()
    """