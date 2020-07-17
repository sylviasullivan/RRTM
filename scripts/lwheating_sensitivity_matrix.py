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

# Matrix to hold the radiative heating profiles and their sensitivities.
H = np.zeros((6,tt.shape[0],lev+1))
dHdT = np.zeros((3,tt.shape[0],lev+1))

# Read in the IWP values.
IWP = np.loadtxt(basedir + '/output/IWP_iterator.txt')
c = 1

for tk in tt:
    if c < 10:
       cstr = '0' + str(c)
    else:
       cstr = str(c)
    lw_Tb  = np.genfromtxt(basedir + '/output/Tb/lwflxatm-Tb_' + str(tk+6) + '.txt')
    lw_Tbp = np.genfromtxt(basedir + '/output/Tb_p/lwflxatm-Tb' + cstr + '__' + str(tk+6) + 'p.txt')
    lw_Tt  = np.genfromtxt(basedir + '/output/Tt/lwflxatm-Tt_' + str(tk) + '.txt')
    lw_Ttp = np.genfromtxt(basedir + '/output/Tt_p/lwflxatm-Tt' + cstr + '__' + str(tk) + 'p.txt')
    lw_IWP  = np.genfromtxt(basedir + '/output/THI/lwflxatm-IWP' + cstr + '_THI.txt')
    lw_IWPp = np.genfromtxt(basedir + '/output/THI_p/lwflxatm-IWP' + cstr + '_THIp.txt')

    # First column is all-sky. Second column is clear-sky.
    # The factor -1 is to define outgoing fluxes as positive.
    lwcld_Tb   = -1.*(lw_Tb[0]-lw_Tb[1])
    lwcld_Tbp  = -1.*(lw_Tbp[0]-lw_Tbp[1])
    lwcld_Tt   = -1.*(lw_Tt[0]-lw_Tt[1])
    lwcld_Ttp  = -1.*(lw_Ttp[0]-lw_Ttp[1])
    lwcld_IWP  = -1.*(lw_IWP[0]-lw_IWP[1])
    lwcld_IWPp = -1.*(lw_IWPp[0]-lw_IWPp[1])

    # Calculate the longwave and shortwave heating rates.
    # Factor of 3600 converts from K s-1 to K day-1.
    H[0,c-1] = g/cp*np.gradient(lwcld_Tb,pp_hl*100.)*3600.
    H[1,c-1] = g/cp*np.gradient(lwcld_Tbp,pp_hl*100.)*3600.
    H[2,c-1] = g/cp*np.gradient(lwcld_Tt,pp_hl*100.)*3600.
    H[3,c-1] = g/cp*np.gradient(lwcld_Ttp,pp_hl*100.)*3600.
    H[4,c-1] = g/cp*np.gradient(lwcld_IWP,pp_hl*100.)*3600.
    H[5,c-1] = g/cp*np.gradient(lwcld_IWPp,pp_hl*100.)*3600.

    # Calculate the sensitivity of the LW heating to the 10% Tb/Tt/IWP perturbation.
    dHdT[0] = (H[1] - H[0])/(0.01*(tk+6))
    dHdT[1] = (H[3] - H[2])/(0.01*tk)
    dHdT[2] = (H[5] - H[4])/(0.01*IWP[c-1])

    c = c + 1


# x-tick labels
xtlbl = [[] for i in np.arange(3)]
xtlbl[0] = tb
xtlbl[1] = tt
xtlbl[2] = [round(i*1000,2) for i in IWP]

for i in np.arange(3):
    print('min / max of sens matrix row ' + str(i) + ': ' + str(np.nanmin(dHdT[i])) + '  ' + str(np.nanmax(dHdT[i])))

fig, ax = plt.subplots(nrows=1,ncols=3,figsize=(14,5))

# Some formatting factors for the plot.. fontsize, tick intervals, and iterator.
fs = 13
c = 0
mm = np.array([[-0.7,0.7],[-0.7,0.7],[-2,2]])
lwsw = ['LW','SW']
xlbl = [r'$T_b$ with $T_t$ = 200 K',r'$T_t$ with $T_b$ = 237 K','IWP [g kg$^{-1}$']
#lbl = [r'K day$^{-1}$ K$^{-1}$',r'K day$^{-1}$ K$^{-1}$',r'K day$^{-1}$ per g kg$^{-1}$']
lbl = ['','','K day$^{-1}$ per g m$^{-2}$']

yt_tt_hl = np.array([int(t) for t in tt_hl[tropopause+3:49:5]])

for j in np.arange(3):
    sns.heatmap(dHdT[j].T,cmap=cm.bwr,ax=ax[j],center=0.,vmin=mm[c,0],vmax=mm[c,1],\
               cbar_kws={'label':lbl[c]})

    # Adjust tick labels
    if j == 0:
       ax[j].set_ylabel('Temperature [K]',fontsize=fs)
    ax[j].set_xlabel(xlbl[j],fontsize=fs)

    ax[j].set_ylim([49,tropopause+3]) # ylim is done based on indices not values.
    ax[j].set_yticks(np.arange(tropopause+3,49,5))
    ax[j].set_yticklabels(yt_tt_hl,rotation=45,fontsize=fs-2)
    ax[j].set_xticks(np.arange(0,len(xtlbl[j]),3))  # Label every 2nd tick.
    ax[j].set_xticklabels(xtlbl[j][::3],rotation=45,fontsize=fs-2)

    # Where is 200 K?
    ax[j].plot([0,len(xtlbl[j])],[t200+tropopause,t200+tropopause],linewidth=1,color='k',linestyle='--')
    # Where is 237 K?
    ax[j].plot([0,len(xtlbl[j])],[t237+tropopause+3,t237+tropopause+3],linewidth=1,color='k',linestyle='--')

    c = c + 1

fig.savefig('../figures/lwheating_sensitivity_matrix.pdf',bbox_inches='tight')
plt.show()
