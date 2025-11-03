import numpy as np
import matplotlib.pyplot as plt
import math
from scipy import interpolate
import seaborn as sns
from mpl_toolkits.axes_grid1 import make_axes_locatable
from astropy import constants as const
import random
import scipy 

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
 
files = ['los2048_n5000_z5.200.dat', 'los2048_n5000_z6.100.dat', 'los2048_n5000_z7.300.dat', 'los2048_n5000_z9.000.dat', 'los2048_n5000_z12.900.dat']
#read in text file of file names
#lines = np.genfromtxt('filenames.txt',dtype='str')
#print(len(lines))
#choose redshift slice from o to 200
#z_file = 156

"""_
 ------------------------------------------------------------------------------------------------------------------------
 Read in binary file data 
 ------------------------------------------------------------------------------------------------------------------------
"""
all_mfp = np.empty([len(files),5000])
all_z = np.array([])
for k in range(0, len(files)):
    # File directory and name for lyman alpha spectra things
    file  = files[k]

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

    #print('This is for redshift' + str(ztime))

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

    tau_Lya = cross_sec*n_HI*delta_r

    for j in range (0, 5000):
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
        #print(n*delta_r)
        #print(n)

    #convert to pMpc from pcm 
    proper_mfp = (mfp*1e-2 /(MPC))
    #proper_mfp = mfp /3.086e+24
    all_mfp[k,:] = proper_mfp
    all_z = np.append(all_z, ztime)

        

    #define points from literature - in pMpc
    obs_mfp = np.array([[6, 5.1, 5.16, 4.86, 4.56],[0.75, 9.09, 10.3, 15.1, 22.2], [0.45, 1.28, 1.6, 1.8, 2.3], [0.65, 1.62, 1.6, 1.8, 2.3]])

   
    
    #histogram has showed mfp is biomodal so I should plot both as seperate points on the histogram (gonna have to think of a way to automate this)
    """
    plt.subplots(1,1)
    plt.errorbar(obs_mfp[0,:], obs_mfp[1,:],yerr=obs_mfp[2:4,:], fmt="o", label = 'Observations')
    #plt.scatter(ztime, np.mean(mfp_1), color ='r')
    #plt.scatter(ztime, np.mean(mfp_2), color ='k')
    plt.errorbar(ztime, np.mean(proper_mfp), yerr=np.std(proper_mfp), fmt="o", label = 'Sherwood Mean')
    plt.errorbar(ztime, np.median(proper_mfp), yerr=np.std(proper_mfp), fmt="o", label = 'Sherwood Median')
    #plt.errorbar(ztime, np.mean(proper_mfp_HI), yerr=np.ptp(proper_mfp_HI), fmt="o", label = 'Sherwood HI')
    plt.yscale('log')
    plt.ylabel(r'$\lambda_{mfp}$ (pMpc)')
    plt.xlabel('z')
    plt.legend()
    #plt.show()
    

    # Plot the transmitted flux along a selected line of sight
    PLOTLOS = 4999 # select line of sight for plotting (0->numlos-1)

    F = np.array([])
    tau_eff = np.array([])
    for i in range(0, int(numlos -1 )):
        tau_eff = np.append(tau_eff, -np.log(np.sum(np.exp(-tau_Lya[PLOTLOS*nbins[0] : (PLOTLOS+1)*nbins[0]]))/((PLOTLOS+1)*nbins[0] - (PLOTLOS)*nbins[0])))
        
        F =np.append(F,np.sum(np.exp(-tau_Lya[PLOTLOS*nbins[0] : (PLOTLOS+1)*nbins[0]]))/((PLOTLOS+1)*nbins[0] - (PLOTLOS)*nbins[0]))
    
    plt.figure(3)
    plt.plot(velaxis,-(np.exp(-tau_Lya[PLOTLOS*nbins[0] : (PLOTLOS+1)*nbins[0]])-1),color='blue')
    plt.xlabel(r'$\rm v_{H}\,[km\, s^{-1}]$',fontsize=12)
    plt.ylabel(r'$\rm F=e^{-\tau}$',fontsize=12)
    plt.title('z = ' + str(ztime))
    #plt.show()

    print('at z = ' + str(ztime) + ' tau_eff = ' +str(np.mean(tau_eff)) + " <F> = " + str(np.mean(F)))
    """ 
   


count, bins_count = np.histogram(all_mfp[0,:], bins=10)

# finding the PDF of the histogram using count values
pdf = count / sum(count)
 
plt.figure()
print(all_mfp[0,:])

for i in range(0,len(all_mfp)):
    #plt.hist(all_mfp[i,:])
    """
    y, binEdges = np.histogram(all_mfp[i,:], bins=100)
    bincenters = 0.5 * (binEdges[1:] + binEdges[:-1])
    plt.plot(bincenters, y, '-', c='black')
    """    
    count, bins_count = np.histogram(all_mfp[i,:], bins=100)
    pdf = count / sum(count)
    plt.plot(bins_count[1:], pdf, label=str(all_z[i]))

plt.xlabel(r'$\lambda_{mfp}$ (pMpc)')
plt.ylabel('PDF')
plt.xscale('log')
plt.legend()
plt.show()
for i in range(0,len(all_mfp)):
    print('at z = ' + str(all_z[i]) + ' <mfp> = ' +str(np.mean(all_mfp[i,:])))