#make first paper overview plot
import pandas as pd
import numpy as np
from plotting_dictonary import dictionary
import matplotlib.pyplot as plt
import matplotlib

#Paper specific Matploblib settings
import matplotlib.font_manager
matplotlib.style.use('paper_params.mplstyle')
import matplotlib as mpl

def log_lims(value, up_err, low_err):
	err_down = 10**(value - low_err) - 10**(value)
	err_up =  10**(value + up_err) - 10**(value)
	return err_down, err_up

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

nrows = 2

ncols = 2

#import data for graph

"""average transmisssion"""
#load in data 
base = '/home/ppxjf3/first_year/Sherwood_Relics/'
F_zr57 = np.loadtxt(base + 'Optical_depth/planck1_40_2048_RTzrfit_los.txt')
F_zr57_T40 = np.loadtxt(base + 'Optical_depth/planck1_160_2048_RTzrfit_T40_los.txt')
F_zr53 = np.loadtxt(base + 'Optical_depth/planck1_40_2048_RTzr53.txt')
F_zr57_low_res = np.loadtxt(base + 'Optical_depth/planck1_40_512_RTzrfit.txt')


#define observational data 
#Bosman+22
z_b22= np.arange(4.8,6.2, 0.1)
F_b22 = [0.194, 0.171, 0.1581, 0.1428, 0.1222, 0.1031, 0.0801, 0.0591, 0.0447, 0.0256, 0.0172, 0.0114, 0.0089, 0.0088, 0.0047]
err_up = [0.018,0.014,0.0082,0.0068,0.0046,0.0056,0.0061,0.0039,0.0033,0.0031,0.0022,0.0029,0.0033,0.0082,0.0045]
err_down = [0.015,0.014,0.0089, 0.0054, 0.0054, 0.005, 0.0048, 0.0035, 0.0036, 0.0029, 0.0028, 0.003, 0.0029, 0.0074, 0.0044]

#bosman+18
z_b_18 = [5.0, 5.2, 5.4, 5.6, 5.8, 6.0]
F_b_18 = [0.135, 0.114, 0.084, 0.05, 0.023, 0.0072]
err_b18 = [0.012, 0.006, 0.005, 0.005, 0.004, 0.0018]

#eilers 18
z_e18 = [4.25, 4.5, 4.75, 5.0, 5.25, 5.5, 5.75, 6.0, 6.25]
F_e18 = [0.3595, 0.2927, 0.1944, 0.1247, 0.0795, 0.0531, 0.0182, 0.0052, -0.0025]
err_e18 = [0.0112, 0.019, 0.015, 0.0132, 0.0078, 0.0058, 0.0045, 0.0043, 0.0007]

#becker 2013
z_b13 = [4.55 , 4.65 , 4.75 ,  4.85 ]
F_b13 = [0.3009, 0.2881, 0.2419, 0.2225]
err_b13 = [ 0.0104,  0.0117, 0.0201, 0.0151]

""" MFP """
#load in data for 5000 itts
mfp_zr57 = np.loadtxt(base + 'MFP_plots/itter105/planck1_40_2048_RTzrfit_los_mfp_N100000.txt')
idx = np.argsort(mfp_zr57[0,:])
mfp_zr57[0,:] = mfp_zr57[0,idx]
mfp_zr57[1,:] = mfp_zr57[1,idx]
mfp_zr57_T40 = np.loadtxt(base + 'MFP_plots/itter105/planck1_160_2048_RTzrfit_T40_los_mfp_N100000.txt')
idx = np.argsort(mfp_zr57_T40[0,:])
mfp_zr57_T40[0,:] = mfp_zr57_T40[0,idx]
mfp_zr57_T40[1,:] = mfp_zr57_T40[1,idx]
mfp_zr53 = np.loadtxt(base + 'MFP_plots/planck1_40_2048_RTzr53_mfp_N10000.txt')
idx = np.argsort(mfp_zr53[0,:])
mfp_zr53[0,:] = mfp_zr53[0,idx]
mfp_zr53[1,:] = mfp_zr53[1,idx]
mfp_zr57_low_res = np.loadtxt(base + 'MFP_plots/planck1_40_512_RTzrfit_mfp_N10000.txt')
idx = np.argsort(mfp_zr57_low_res[0,:])
mfp_zr57_low_res[0,:] = mfp_zr57_low_res[0,idx]
mfp_zr57_low_res[1,:] = mfp_zr57_low_res[1,idx]

#mfp data Gaikwad
h100 = 0.678
x_g23 = np.array([4.9, 5.0, 5.1, 5.2, 5.3, 5.4, 5.5, 5.6, 5.7, 5.8, 5.9, 6.0])
y_g23 = np.array([50.119, 57.412, 46.452, 40.832, 34.041, 29.242, 28.907, 22.961, 16.596, 13.183, 10.471, 8.318])/((1+x_g23)*h100)
yerr_g23 = np.array([[19.919, 33.058], [21.104, 38.088], [18.909, 31.173], [15.713,22.264], [12.163, 18.44], [10.187, 15.427], [11.124, 10.904], [7.826, 11.712], [7.476, 9.707], [5.769, 9.205], [4.976, 9.027], [4.052, 7.531]]).T/((1+x_g23)*h100)
xerr_g23 = 0.05

#Satyavolu 23
x_s23 = np.array([6.00])
mfp_s23 = np.array([0.90])
mfp_err_s23 = np.array([[0.40,0.66]]).T

#Zhu 23
x_z23_mfp = np.array([5.08, 5.31, 5.65, 5.93])
mfp_z23 = np.array([9.33, 5.40, 3.31, 0.81]) #pMpc
mfp_err_z23 = np.array([[1.80,2.06],[1.40,1.47],[1.34,2.74],[0.48,0.73]]).T

#mfp sim data
mfp_keat20 = np.loadtxt(base + 'other_sim_mfp_data/z_mfp_lowcmbtau.txt')
mfp_lew22 = np.loadtxt(base + 'other_sim_mfp_data/lewis_22.txt')
mfp_gara21 = pd.read_pickle(base + 'other_sim_mfp_data/mfp_Thesan-1.pkl')  


"""Gamma"""
z_gamma = np.loadtxt(base + 'gamma_data/gamma_redshifts.txt')
av_gamma1 = np.loadtxt(base + 'gamma_data/zr57_xh1.txt')
av_gamma2 = np.loadtxt(base + 'gamma_data/zr57_T40_xh1.txt')
av_gamma3 = np.loadtxt(base + 'gamma_data/zr53_xh1.txt')
av_gamma4 = np.loadtxt(base + 'gamma_data/zr57_lowres_xh1.txt')

#obvs data 
#Gaikward 2023
xg_g23 = [4.9, 5.0, 5.1, 5.2, 5.3, 5.4, 5.5, 5.6, 5.7, 5.8, 5.9, 6.0]
yg_g23 = [0.501, 0.557, 0.508, 0.502, 0.404, 0.372, 0.344, 0.319, 0.224, 0.178, 0.151, 0.145]
xerrg_g23 = 0.05
yerrg_g23 = np.array([(0.232, 0.275), (0.218, 0.376), (0.192, 0.324), (0.193, 0.292), (0.147, 0.272), (0.126, 0.217), (0.13, 0.219), (0.12, 0.194), (0.112, 0.223), (0.078, 0.194),(0.079, 0.151), (0.087, 0.157)]).T

#Calverley 2011
yg_c11 = [10**-0.15, 10**-0.84] #[0.08, 0.08, 0.03, 0.10, 0.012, 0.23, 0.43, 0.14, 0.2, 1.82, 0.48, 0.95, 1.02, 3.89, 1.88]
xg_c11 = [5.04, 6.09] #[6.4189, 6.308, 6.247, 6.02, 6.016, 5.82, 5.81, 5.41, 5.33, 5.2, 5.09, 4.967, 5.886, 4.876, 4.588]
value = np.array([-0.15, -0.84])
up_err = np.array([0.16, 0.18])
down_err = np.array([0.16, 0.18])

yerrgdown_c11, yerrgup_c11 = log_lims(value, up_err, down_err)
yerrg_c11 = np.array([-yerrgdown_c11, yerrgup_c11])
print(yerrg_c11)
xerrg_c11 = np.array([(0.44,0.36), (0.29,0.31)])

#Becker 2013
xg_b13 = [4.0, 4.4, 4.75]
value = np.array([-0.072, -0.019, -0.029])
yg_b13 = [10**-0.072, 10**-0.019, 10**-0.029]
up_err = np.array([0.135, 0.140, 0.156])
down_err = np.array([0.117, 0.122, 0.147])

yerrgdown_b13, yerrgup_b13 = log_lims(value, up_err, down_err)
yerrgdown_b13 = - yerrgdown_b13

"""XHI"""
zxh1 = np.loadtxt(base + 'xh1_data/redshifts.txt')
idx = np.argsort(zxh1)
zxh1 = zxh1[idx]
av_H1frac1 = np.loadtxt(base + 'xh1_data/zr57_xh1.txt')[idx]
av_H1frac2 = np.loadtxt(base + 'xh1_data/zr57_T40_xh1.txt')[idx]
av_H1frac3 = np.loadtxt(base + 'xh1_data/zr53_xh1.txt')[idx]
av_H1frac4 = np.loadtxt(base + 'xh1_data/zr57_lowres_xh1.txt')[idx]
print(zxh1)
#obvs data 
#McGreer 2015
y_mg15 = np.array([0.06, 0.04, 0.58])
x_mg15 = np.array([5.9, 5.6, 6.1])
yerr_mg15 = np.array([[0.05, 0.05, 0.20], [0.0, 0.0, 0.0]])

#Greig 2018
y_g18 = 0.4
x_g18 = 7.1
yerr_g18 = np.array([(0.19,0.21)]).T

#Greig 2022
x_g22 = np.array([7.29])
y_g22 = np.array([0.49])
yerr_g22 = np.array([(0.11,0.11)]).T

#Jin 2023
x_j23 = np.array([5.5, 5.69, 5.89, 6.1, 6.3, 6.49, 6.69])
y_j23 = np.array([0.076, 0.16, 0.28, 0.68, 0.78, 0.86, 0.94])
yerr_j23 =  np.array([(0.08,0.00), (0.13,0.0), (0.07, 0.0), (0.055,0.0), (0.05,0.0), (0.03,0.0), (0.06,0.0)]).T

#Zhu 2022
y_z22 = np.array([0.167, 0.281,0.045])
x_z22 = np.array([5.76, 5.94, 5.53])
err_z22 = np.array([(0.09,0.0), (0.09,0.0), (0.03,0.0)]).T

#Umeda 2023
x_u23 = np.array([7.12,7.44,8.28,10.28])
y_u23 = np.array([0.54,0.69,0.92,0.94])
yerr_u23 = np.array([(0.13,0.54),(0.30,0.38),(0.08,0.56),(0.06,0.41)]).T
xerr_u23 = np.array([(0.08, 0.06),(0.24,0.34),(0.44,0.41),(1.40,1.21)]).T

""" --------- Emissivity ---------"""

basee = '/home/ppxjf3/Downloads/'
#MFP data
data = np.loadtxt(basee + 'ndot_160_2048.txt')
ndot_z_160 = data[:,0]
ndot_160 = data[:,1]

data = np.loadtxt(basee + 'ndot_L40N2048_bestfit.txt')
ndot_z_40 = data[:,0]
ndot_40 = data[:,1]*(10**-50)

data = np.loadtxt(basee + 'z_em_lowcmbtau.txt')
ndot_z_k20 = data[:,0]
ndot_k20 = data[:,1]*(10**-50)
#make graph 
# axs is a numpy array with dimension (nrows, ncols)

fig, ax = plt.subplots(figsize=(fig_width*1.5, fig_height*1.5))

#av transmission      
ax5=fig.add_subplot(323)                 
obs1 = ax5.errorbar(z_b22, F_b22, yerr = np.array([err_down, err_up]),fmt=dictionary['obs1']['fmt'], label = 'Bosman+22', capsize=6, color =dictionary['obs1']['c'])
obs2 = ax5.errorbar(z_b_18, F_b_18, yerr = err_b18, fmt=dictionary['obs2']['fmt'], label = 'Bosman+18', capsize=6, color =dictionary['obs2']['c'])
obs3 = ax5.errorbar(z_e18, F_e18, yerr = err_e18, fmt=dictionary['obs3']['fmt'], label = 'Eilers+18', capsize=6, color =dictionary['obs3']['c'])
#obs4 = ax[0,0].errorbar(z_b13, F_b13, yerr = err_b13, fmt=dictionary['obs4']['fmt'], label = 'Becker+13', capsize=6, color =dictionary['obs4']['c'])

sim_1, = ax5.plot(F_zr57[0,:], F_zr57[1,:], linestyle = dictionary['zr57']['ls'] , color = dictionary['zr57']['c'], linewidth=1.5)
sim_2, = ax5.plot(F_zr57_T40[0,:], F_zr57_T40[1,:], linestyle = dictionary['zr57_T40']['ls'] , color = dictionary['zr57_T40']['c'])

ax5.set_ylabel(r'$\langle F \rangle $')
ax5.set_xlabel('z')
ax5.legend(frameon=False)
#first_legend = ax5.legend(handles=[obs1, obs3, obs2], loc='upper right', frameon=False)
#ax5.add_artist(first_legend)
#ax5.legend(handles=[sim_1, sim_2], loc='lower left', frameon=False)
ax5.set_xlim(left = 4.7, right = 6.5)
ax5.set_ylim(top = 0.25)

#emissivity
arrow_size =0.1
ax4=fig.add_subplot(324)
ax4.plot(ndot_z_160, ndot_160, linestyle = dictionary['zr57_T40']['ls'] , color = dictionary['zr57_T40']['c'])
ax4.plot(ndot_z_40, ndot_40, linestyle = dictionary['zr57']['ls'] , color = dictionary['zr57']['c'], linewidth=1.5)
#ax4.plot(ndot_z_k20, np.log10(ndot_k20), linestyle = '-.', color = 'brown')
ax4.set_ylabel(r'$\dot{n}/10^{50} \rm s^{-1} cMpc^{-3}$')
ax4.set_xlabel('z')
ax4.legend(frameon = False)

ax4.set_xlim(left=4.9, right=12)
ax4.set_ylim(bottom=0.80, top =12)
ax4.set_yscale('log')


#xHI

arrow_size =0.1
ax3=fig.add_subplot(325)
sim1, = ax3.plot(zxh1, av_H1frac1, linestyle =  dictionary['zr57']['ls'], color =  dictionary['zr57']['c'], linewidth=1.5)
sim2, = ax3.plot(zxh1, av_H1frac2, linestyle = dictionary['zr57_T40']['ls'], color = dictionary['zr57_T40']['c'])
#sim3, = ax[1,0].plot(zxh1, av_H1frac2, label= dictionary['zr53']['label'], linestyle= dictionary['zr53']['ls'] , color = dictionary['zr53']['c'])
#sim4, = ax[1,0].plot(zxh1, av_H1frac4, label= dictionary['zr57_low_res']['label'], linestyle= dictionary['zr57_low_res']['ls'] , color = dictionary['zr57_low_res']['c'])

obs5 = ax3.errorbar(x_g18, 1-y_g18, yerr= yerr_g18, linestyle ='None', label = 'Greig+17', capsize=6, color = dictionary['obs5']['c'],  fmt=dictionary['obs5']['fmt'])
obs6 = ax3.errorbar(x_mg15, 1-y_mg15, yerr=arrow_size, linestyle ='None', label = 'McGreer+15',   lolims=True, color = 'b', fmt="X")
obs6pl, obs6cap, obs6bar = ax3.errorbar(x_mg15, 1-y_mg15, yerr= yerr_mg15, linestyle ='None', label = 'McGreer+15',   uplims=True, color = 'b', fmt="X")
for cap in obs6cap:
    cap.set_marker('_')
obs2 = ax3.errorbar(x_j23, 1- y_j23, linestyle ='None', yerr = arrow_size, label = 'Jin+23', lolims=True, color = dictionary['obs2']['c'], fmt=dictionary['obs2']['fmt'])
obs2pl, obs2cap, obs2bar = ax3.errorbar(x_j23, 1- y_j23, linestyle ='None', yerr = yerr_j23, label = 'Jin+23', uplims=True, color = dictionary['obs2']['c'], fmt=dictionary['obs2']['fmt'])
for cap in obs2cap:
    cap.set_marker('_')
obs4 = ax3.errorbar(x_g22, 1- y_g22, yerr= yerr_g22, linestyle ='None', capsize=6, label = 'Greig+22', color = dictionary['obs4']['c'],  fmt=dictionary['obs4']['fmt'])
obs3 = ax3.errorbar(x_z22, 1- y_z22, linestyle ='None', yerr = arrow_size, label = 'Zhu+22', lolims=True, color = dictionary['obs3']['c'], fmt=dictionary['obs3']['fmt'])
obs3pl, obs3cap, obs3bar = ax3.errorbar(x_z22, 1- y_z22, linestyle ='None', yerr = err_z22, label = 'Zhu+22', uplims=True, color = dictionary['obs3']['c'], fmt=dictionary['obs3']['fmt'])
for cap in obs3cap:
    cap.set_marker('_')
obs1 = ax3.errorbar(x_u23, 1- y_u23, linestyle ='None', yerr = yerr_u23, xerr=xerr_u23, label = 'Umeda+23', capsize=6, color = dictionary['obs1']['c'], fmt=dictionary['obs1']['fmt'])
ax3.set_ylabel(r'1 - $\langle x_{\mathrm{HI}} \rangle$')
ax3.set_xlabel('z')
first_legend = ax3.legend(handles=[obs1, obs2, obs3, obs4, obs5, obs6], loc='upper right', frameon=False)
ax3.add_artist(first_legend)
#ax[1,0].legend(handles=[sim_1, sim_2, sim4], loc='lower left', frameon=False)
ax3.set_xlim(left = 4.7, right=12)

#Gamma
ax2=fig.add_subplot(326)
sim1, = ax2.plot(z_gamma, av_gamma1, linestyle =  dictionary['zr57']['ls'], color =  dictionary['zr57']['c'], linewidth=1.5)
sim2, = ax2.plot(z_gamma, av_gamma3, linestyle = dictionary['zr57_T40']['ls'], color = dictionary['zr57_T40']['c'])
#sim3, = ax[1,1].plot(z_gamma, av_gamma2, label= dictionary['zr53']['label'], linestyle= dictionary['zr53']['ls'] , color = dictionary['zr53']['c'])
#sim4, = ax[1,1].plot(z_gamma, av_gamma4, label= dictionary['zr57_low_res']['label'], linestyle= dictionary['zr57_low_res']['ls'] , color = dictionary['zr57_low_res']['c'])

obs1 = ax2.errorbar(xg_g23, yg_g23, yerr= yerrg_g23, xerr=xerrg_g23, linestyle ='None', label = 'Gaikwad+23', capsize=6, color = dictionary['obs1']['c'],  fmt= dictionary['obs1']['fmt'])
obs2 = ax2.errorbar(xg_b13, yg_b13, yerr= np.array([yerrgdown_b13, yerrgup_b13]),  linestyle ='None', label = 'Becker+13', capsize=6, color = dictionary['obs2']['c'],  fmt= dictionary['obs2']['fmt'])
obs3 = ax2.errorbar(xg_c11, yg_c11, yerr= yerrg_c11, xerr = xerrg_c11, linestyle ='None', label = 'Calverley+11', capsize=6, color = dictionary['obs3']['c'],  fmt= dictionary['obs3']['fmt'])

ax2.set_ylabel(r'$\langle \Gamma_{\mathrm{HI}} \rangle / 10^{-12} \rm s^{-1}$' )
ax2.set_xlabel('z')
first_legend = ax2.legend(handles=[obs1, obs2, obs3], loc='upper right', frameon=False)
ax2.add_artist(first_legend)
#ax[1,1].legend(handles=[sim1, sim2, sim4], loc='lower left', frameon=False)
ax2.set_xlim(left = 4.7, right=7)
ax2.set_ylim(bottom=0.05)
ax2.set_yscale('log')


#mfp
ax1=fig.add_subplot(311)
obs_4 = ax1.errorbar(np.array([4.56, 4.86, 5.16]),np.array([22.2, 15.1, 10.3]), yerr = np.array([2.3, 1.8, 1.6]), fmt=dictionary['obs4']['fmt'] , label = 'Worseck+14', capsize=6, color = dictionary['obs4']['c'])
obs_3 = ax1.errorbar(np.array([6.0, 5.1]),np.array([0.75, 9.09]), yerr = np.array([[0.45, 1.28], [0.65,  1.62]]), fmt=dictionary['obs3']['fmt'], label = 'Becker+21', capsize=6, color = dictionary['obs3']['c']) #lower err first 
obs_5 = ax1.errorbar(np.array([4.07, 4.22]),np.array([33.0, 28.1]), yerr = np.array([3.5, 2.9]), fmt=dictionary['obs5']['fmt'], label = 'Prochaska+09', capsize=6, color = dictionary['obs5']['c']) #lower err first 
obs_2 = ax1.errorbar(x_g23, y_g23, yerr = yerr_g23, xerr = xerr_g23, fmt=dictionary['obs2']['fmt'], label = 'Gaikwad+23', capsize=6, color =dictionary['obs2']['c'])
obs_1 = ax1.errorbar(x_z23_mfp, mfp_z23, yerr = mfp_err_z23, label = 'Zhu+23', capsize=6, fmt=dictionary['obs1']['fmt'], color=dictionary['obs1']['c'])   
obs_5 = ax1.errorbar(x_s23, mfp_s23, yerr = mfp_err_s23, label = 'Satyavolu+23', capsize=6, fmt='*', color='darkmagenta')  

sim_1, = ax1.plot(mfp_zr57[0,:], mfp_zr57[1,:], linestyle = dictionary['zr57']['ls'] , color = dictionary['zr57']['c'], linewidth=1.5)
sim_2, = ax1.plot(mfp_zr57_T40[0,:], mfp_zr57_T40[1,:], linestyle = dictionary['zr57_T40']['ls'] , color = dictionary['zr57_T40']['c'])
#sim_3, = ax[0,1].plot(mfp_zr53[0,:], mfp_zr53[1,:], label= 'zr53', linestyle= dictionary['zr53']['ls'], color =  dictionary['zr53']['c'])
#sim_4 ,= ax[0,1].plot(mfp_zr57_low_res[0,:], mfp_zr57_low_res[1,:], label= 'zr57_low_res', linestyle= dictionary['zr57_low_res']['ls'], color =  dictionary['zr57_low_res']['c'])
sim_sher = ax1.plot([],[], label= dictionary['zr57']['label'] , linestyle = dictionary['zr57']['ls'] , color = dictionary['zr57']['c'], linewidth=1.5)
sim_sher2 = ax1.plot([],[], label= dictionary['zr57_T40']['label'] , linestyle = dictionary['zr57_T40']['ls'] , color = dictionary['zr57_T40']['c'])
sim_k20 = ax1.plot(mfp_keat20[:,0], mfp_keat20[:,1]/h100, linestyle = '-.', color = 'brown', label='Keating+20b')
sim_l22 = ax1.plot(mfp_lew22[0,:], mfp_lew22[1,:], linestyle = '--', color = 'blue', label='Lewis+22')
sim_g21 = ax1.plot(mfp_gara21['redshift'].values, mfp_gara21['mfp'].values, linestyle = 'solid', color = 'gray', label='Garaldi+22')

curves = sim_l22 + sim_g21 + sim_k20 + sim_sher + sim_sher2
ax1.set_yscale('log')
ax1.set_ylabel(r'$\lambda_{\mathrm{mfp}}$ [pMpc]')
ax1.set_xlabel('z')
first_legend = ax1.legend(handles=[obs_1, obs_5, obs_2, obs_3, obs_4], loc='lower left', frameon=False)
ax1.add_artist(first_legend)
labs = [curve.get_label() for curve in curves]
second_legend = ax1.legend(curves, labs, loc='upper right', frameon=False)
ax1.set_xlim(left = 4.7, right = 7)
ax1.set_ylim(bottom = 0.2, top=30)
fig.delaxes(ax)

plt.savefig(base+'figure1_new.pdf', bbox_inches='tight') #, pad_inches=0.02)  # remove whitespace
plt.show()

nrows = 1

ncols = 1

data = np.loadtxt(basee + 'z_meanflux_lowcmbtau.txt')
F_z_k20 = data[:,0]
F_k20 = data[:,1]
fig, ax = plt.subplots(figsize=(fig_width, fig_height),
                       nrows=nrows, ncols=ncols)

obs1 = ax.errorbar(z_b22, F_b22, yerr = np.array([err_down, err_up]),fmt=dictionary['obs1']['fmt'], label = 'Bosman+22', capsize=6, color =dictionary['obs1']['c'])
obs2 = ax.errorbar(z_b_18, F_b_18, yerr = err_b18, fmt=dictionary['obs2']['fmt'], label = 'Bosman+18', capsize=6, color =dictionary['obs2']['c'])
obs3 = ax.errorbar(z_e18, F_e18, yerr = err_e18, fmt=dictionary['obs3']['fmt'], label = 'Eilers+18', capsize=6, color =dictionary['obs3']['c'])
#obs4 = ax[0,0].errorbar(z_b13, F_b13, yerr = err_b13, fmt=dictionary['obs4']['fmt'], label = 'Becker+13', capsize=6, color =dictionary['obs4']['c'])

sim_1, = ax.plot(F_zr57[0,:], F_zr57[1,:], label='40-2048', linestyle = dictionary['zr57']['ls'] , color = dictionary['zr57']['c'], linewidth=1.5)
sim_2, = ax.plot(F_zr57_T40[0,:], F_zr57_T40[1,:], label='160-2048', linestyle = dictionary['zr57_T40']['ls'] , color = dictionary['zr57_T40']['c'])
ax.plot(F_z_k20, F_k20, color='brown', linestyle='dashdot', label='Keating+20b')
ax.set_ylabel(r'$\langle F \rangle $')
ax.set_xlabel('z')
ax.legend(frameon=False)
#first_legend = ax5.legend(handles=[obs1, obs3, obs2], loc='upper right', frameon=False)
#ax5.add_artist(first_legend)
#ax5.legend(handles=[sim_1, sim_2], loc='lower left', frameon=False)
ax5.set_xlim(left = 4.7, right = 6.5)
ax5.set_ylim(top = 0.25)

plt.tight_layout()
plt.savefig(base + 'ave_flux.pdf', dpi=330, bbox_inches='tight')
plt.show()
