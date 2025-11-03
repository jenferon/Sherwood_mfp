import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
from astropy import constants as const

#global constants
MPC        = 3.08568025e+22  # m
H0         = 1.0e5/MPC       # s-1
kPC        = const.kpc.value #m
GRAVITY    = const.G.value   #m3 / (kg s2)
PROTONMASS = const.m_p.value #kg
kb = const.k_B.value # J / (K)

def number_density_H(h100,overdensity, omegab, XHfrac, redshift):
    #Get nH(z) at the critical density in proper cgs units
    
    rhoc       = 3.0 * (H0*h100)**2/(8.0 * np.pi * GRAVITY) #kg m-3
    
    n_H     = rhoc*omegab*(overdensity)*XHfrac*(1.0 + redshift)**3.0 / PROTONMASS #m-3
    
    return n_H

def find_nearest(array, value):
		array = np.asarray(array)
		idx = (np.abs(array - value)).argmin()
		return idx  


def pdf_calc (bin_min, bin_max, nbins, data, dr):
	step = (bin_max - bin_min)/nbins
	#print(step)
	bins = np.arange(bin_min, bin_max, step)
	N = np.array([])
	for j in range (0, nbins-1):
		
		Ni =np.where((data>= bins[j]) & (data < bins[j+1]))
		N = np.append(N,len(Ni[0]))
	
	pdf = N/(step*dr)
	centers=0.5*(bins[:-1]+bins[1:])
	return pdf, bins, centers


def normalisation_factor(v, itter, h100):
    R = v/(H0*h100) #km
    R_m = R*100 #m
    R_MPc = R_m/MPC
    return R_MPc*itter

def trap_integ(x, y):
    if len(x) == 0:
        return 0
    idx = np.argsort(x)
    x = x[idx]
    y = y[idx]
    sum = 0 
    for i in range(1,len(y)-2):
        sum += ((y[i-1] + y[i])/2 )*(x[i]- x[i-1])
    return sum

def logmean(x):
    return np.log10(np.mean(x))

def jeans_scale(temp, z):
    Y=0.24
    z = float(z)
    
    ne_nH = 1 +  (Y / (4*(1-Y)))
    mu = ((1-Y)*ne_nH+ Y/4)**-1
    T = np.mean(temp)/(10**2)
    lambda_j = 9.8*np.sqrt(10*(1.22/mu)*T*(10/(1+z)))
    return lambda_j


def column_density(window=0.5*1000, window_func='proper', thresh= 0.5, z = ['6.0', '10.0', '8.0', '7.0', '5.4'], folder = 'planck1_40_2048_RTzrfit_los'):

    N_HI_j = np.array([])
    N_HI_c = np.array([])
    N_HI_p = np.array([])
    for jj in z:
        file = '/home/ppxjf3/' + folder + '/los2048_n5000_z' + str(jj)+'00.dat'
        """_
            ------------------------------------------------------------------------------------------------------------------------
            Read in binary file data 
            ------------------------------------------------------------------------------------------------------------------------
        """
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
        --------------------------------------------------------------------------------------------------------------------------------------------
        Column Density calculation
        -----------------------------------------------------------------------------------------------------------------------------------------------
        """
        itter = 100000
        #define n and nHI and delta r
        n = number_density_H(h100,density, omegab, Xh, ztime)*1e-6 #cm-3
        n_HI = n*H1frac #cm-3
        cm_pos = (posaxis*kPC*1e2)/((1+ztime)*h100)#proper cm
        if window_func == 'proper':
            idx = find_nearest(posaxis/(1+ztime), window)
        if window_func == 'comove':
            window_co = ((1+6.0)/(2))*1000 #cMPc/h
            idx = find_nearest(posaxis, window_co)
        if window_func == 'jeans':
            window_j = jeans_scale(temp_1, ztime)
            idx = find_nearest(posaxis, window_j)
        print(idx)
        delta_r = cm_pos[3] - cm_pos[2]

            #random number seeded so the same random numbers are generated each time
        np.random.seed(50) 
        start_idx = np.random.randint(len(n_HI)-1, size = itter)

                
            
        tile_HI = np.tile(n_HI,2) #incase we sample too close to the end of the array it just cycles round
        tile_pos = np.tile(cm_pos,2*numlos[0])
        tile_Temp = np.tile(temp_1,2)
        tile_density = np.tile(density,2)
        log_NHI = np.array([])


        logHNI_nothresh = np.array([])
        logHNI_thresh = np.array([])
        av_T = np.array([])
        av_density = np.array([])
        log_av_T = np.array([])
        log_av_density = np.array([])

        for j in range(0,len(start_idx)):
            nHI_slice = tile_HI[start_idx[j]:start_idx[j]+idx]
            xHI_slice = H1frac[start_idx[j]:start_idx[j]+idx]
            pos_slice = tile_pos[start_idx[j]:start_idx[j]+idx]
            av_T = np.append(av_T,np.mean(tile_Temp[start_idx[j]:start_idx[j]+idx]))
            av_density = np.append(av_density,np.mean(tile_density[start_idx[j]:start_idx[j]+idx]))
            log_av_T = np.append(log_av_T,np.mean(np.log10(tile_Temp[start_idx[j]:start_idx[j]+idx])))
            log_av_density = np.append(log_av_density,np.mean(np.log10(tile_density[start_idx[j]:start_idx[j]+idx])))


            #empty array for column desnity
            NHI = integrate.trapezoid(nHI_slice, dx=delta_r) #cm^-2 trap_integ(pos_slice,nHI_slice) #
            logHNI_nothresh = np.append(logHNI_nothresh, np.log10(NHI))
            if np.sum(xHI_slice) < thresh:
                
                logHNI_thresh = np.append(logHNI_thresh, np.log10(NHI))
            else:
                continue
            #N_HI = np.append(N_HI, np.log10(NHI))
        if jj == '10.0':
            NHI_100 = logHNI_nothresh
            NHI_100_thresh = logHNI_thresh
            av_xHI_100 = np.mean(H1frac)
            dr_100 = normalisation_factor(velaxis[idx], itter, h100) #MPc
        if jj == '8.0':
            NHI_80 = logHNI_nothresh
            NHI_80_thresh = logHNI_thresh
            av_xHI_80 = np.mean(H1frac)
            dr_80 = normalisation_factor(velaxis[idx], itter, h100) #MPc
        if jj == '7.0':
            NHI_70 = logHNI_nothresh
            NHI_70_thresh = logHNI_thresh
            av_xHI_70 = np.mean(H1frac)
            dr_70 = normalisation_factor(velaxis[idx], itter, h100) #MPc
        if jj == '6.0':
            NHI_60 = logHNI_nothresh
            NHI_60_thresh = logHNI_thresh
            av_xHI_60 = np.mean(H1frac)
            dr_60 = normalisation_factor(velaxis[idx], itter, h100) #MPc
        if jj == '5.4':
            NHI_54 = logHNI_nothresh
            NHI_54_thresh = logHNI_thresh
            av_xHI_54 = np.mean(H1frac)
            dr_54 = normalisation_factor(velaxis[idx], itter, h100) #MPc
            
    return NHI_100, NHI_80, NHI_70, NHI_60, NHI_54, dr_100, dr_80, dr_70, dr_60, dr_54, NHI_100_thresh, NHI_80_thresh, NHI_70_thresh, NHI_60_thresh, NHI_54_thresh

NHI_100, NHI_80, NHI_70, NHI_60, NHI_54, dr_100, dr_80, dr_70, dr_60, dr_54, NHI_100_thresh, NHI_80_thresh, NHI_70_thresh, NHI_60_thresh, NHI_54_thresh = column_density()

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
