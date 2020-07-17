import matplotlib.pyplot as plt
import matplotlib.colors as colors
import os, sys, time
from os.path import dirname
import numpy as np
import seaborn as sns
from matplotlib import cm
exec(open('MidpointNormalize.py').read())

basedir = dirname(os.getcwd())

# Number of levels over which the RRTM calculations are done
lev = 81
# Heat capacity [J kg-1 K-1]
cp = 1.08*10**(3)
# Gravity [m s-2]
g = 9.8
# Cloud top and bottom temperatures.
tt = np.arange(200,232)
tb = np.arange(206,238)

# Pressure at half levels
ellingson = np.genfromtxt(basedir + '/output/tropical_profile_ellingson_250m_formatted_top2bottom.txt')
pp_fl = ellingson[:,1]
pp_hl = np.zeros((82,))
for i in np.arange(1,81):
    pp_hl[i] = (pp_fl[i-1] + pp_fl[i])/2.0
pp_hl[0] = pp_fl[0] - (pp_hl[1] - pp_fl[0])
pp_hl[81] = pp_fl[80] + (pp_fl[80] - pp_hl[80])

# Temperature at half levels
tt_hl = ellingson[:,2]
melting_layer = np.argmin(np.abs(tt_hl - 273.15))
tropopause = np.argmin(tt_hl)
t200 = np.argmin(np.abs(tt_hl[tropopause:]-200))
t237 = np.argmin(np.abs(tt_hl[tropopause:]-237))

fs = 13

# Matrix to hold the radiative heating profiles
H = np.zeros((3,tt.shape[0],lev+1))
c = 1
for tk in tt:
    lw_Tb  = np.genfromtxt(basedir + '/output/Tb/lwflxatm-Tb_' + str(tk+6) + '.txt')
    lw_Tt  = np.genfromtxt(basedir + '/output/Tt/lwflxatm-Tt_' + str(tk) + '.txt')
    if c < 10:
       cstr = '0' + str(c)
    else:
       cstr = str(c)
    lw_IWP = np.genfromtxt(basedir + '/output/THI/lwflxatm-IWP' + cstr + '_THI.txt')

    # First column is all-sky. Second column is clear-sky.
    # The factor -1 is to define outgoing fluxes as positive.
    lwcld_Tb  = -1.*(lw_Tb[0]-lw_Tb[1])
    lwcld_Tt  = -1.*(lw_Tt[0]-lw_Tt[1])
    lwcld_IWP = -1.*(lw_IWP[0]-lw_IWP[1])

    # Calculate the longwave and shortwave heating rates.
    # Factor of 3600 converts from K s-1 to K day-1.
    H[0,c-1] = g/cp*np.gradient(lwcld_Tb,pp_hl*100.)*3600.
    H[1,c-1] = g/cp*np.gradient(lwcld_Tt,pp_hl*100.)*3600.
    H[2,c-1] = g/cp*np.gradient(lwcld_IWP,pp_hl*100.)*3600.

    c = c + 1

# Set the xtick labels to Tb, Tt, and IWP values.
xtlbl = [[] for i in np.arange(3)]
xtlbl[0] = tb
xtlbl[1] = tt
# Factor of 1000 to convert kg kg-1 to g kg-1.
IWP = np.loadtxt(basedir + '/output/IWP_iterator.txt')
xtlbl[2] = [round(i*1000,2) for i in IWP]

fig, ax = plt.subplots(nrows=1,ncols=1,figsize=(5,5))

# Some formatting factors for the plot.. fontsize, tick intervals, and iterator.
fs = 13
xlbl = [r'$T_b$ with $T_t$ = 200 K',r'$T_t$ with $T_b$ = 237 K','IWP [g kg$^{-1}$]']
lbl = [r'K day$^{-1}$','',r'K day$^{-1}$']
xtinterval = 3 # Label every 3rd xtick.
ytinterval = 5 # Label every 5th ytick.
yt_tt_hl = np.array([round(t) for t in tt_hl[tropopause+3:49:ytinterval]])


for j in np.arange(3):
    sns.heatmap(H[j].T,cmap=cm.bwr,ax=ax,center=0.,vmin=-1.,vmax=1.8,cbar_kws={'label':lbl[j]})

    # Adjust tick labels
    if j == 0:
       ax.set_ylabel('Temperature [K]',fontsize=fs)
    ax.set_xlabel(xlbl[j],fontsize=fs)

    ax.set_ylim([49,tropopause+3])                   # ylim is done based on indices not values.
    ax.set_yticks(np.arange(tropopause+3,49,ytinterval))
    ax.set_yticklabels(yt_tt_hl,rotation=45,fontsize=fs-2)
    ax.set_xticks(np.arange(0,H.shape[1],xtinterval))
    ax.set_xticklabels(xtlbl[j][::xtinterval],rotation=45,fontsize=fs-2)

    # Where is 200 K?
    #ax[j].plot([0,len(xtlbl[j])],[t200+tropopause,t200+tropopause],linewidth=1,color='k',linestyle='--')
    # Where is 237 K?
    #ax[j].plot([0,len(xtlbl[j])],[t237+tropopause+3,t237+tropopause+3],linewidth=1,color='k',linestyle='--')

fig.savefig('../figures/lwheating_matrix.pdf',bbox_inches='tight')
plt.show()
