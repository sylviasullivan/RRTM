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
H = np.zeros((6,tt.shape[0],lev+1))
c = 1
for tk in tt:
    lw_Tb  = np.genfromtxt(basedir + '/output/Tb/lwflxatm-Tb_' + str(tk+6) + '.txt')
    sw_Tb  = np.genfromtxt(basedir + '/output/Tb/swflxatm-Tb_' + str(tk+6) + '.txt')
    lw_Tt  = np.genfromtxt(basedir + '/output/Tt/lwflxatm-Tt_' + str(tk) + '.txt')
    sw_Tt  = np.genfromtxt(basedir + '/output/Tt/swflxatm-Tt_' + str(tk) + '.txt')
    if c < 10:
       cstr = '0' + str(c)
    else:
       cstr = str(c)
    lw_IWP = np.genfromtxt(basedir + '/output/THI/lwflxatm-IWP' + cstr + '_THI.txt')
    sw_IWP = np.genfromtxt(basedir + '/output/THI/swflxatm-IWP' + cstr + '_THI.txt')

    # First column is all-sky. Second column is clear-sky.
    # The factor -1 is to define outgoing fluxes as positive.
    lwcld_Tb  = -1.*(lw_Tb[0]-lw_Tb[1])
    swcld_Tb  = -1.*(sw_Tb[0]-sw_Tb[1])
    lwcld_Tt  = -1.*(lw_Tt[0]-lw_Tt[1])
    swcld_Tt  = -1.*(sw_Tt[0]-sw_Tt[1])
    lwcld_IWP = -1.*(lw_IWP[0]-lw_IWP[1])
    swcld_IWP = -1.*(sw_IWP[0]-sw_IWP[1])

    # Calculate the longwave and shortwave heating rates.
    # Factor of 3600 converts from K s-1 to K day-1.
    H[0,c-1] = g/cp*np.gradient(lwcld_Tb,pp_hl*100.)*3600.
    H[1,c-1] = g/cp*np.gradient(swcld_Tb,pp_hl*100.)*3600.
    H[2,c-1] = g/cp*np.gradient(lwcld_Tt,pp_hl*100.)*3600.
    H[3,c-1] = g/cp*np.gradient(swcld_Tt,pp_hl*100.)*3600.
    H[4,c-1] = g/cp*np.gradient(lwcld_IWP,pp_hl*100.)*3600.
    H[5,c-1] = g/cp*np.gradient(swcld_IWP,pp_hl*100.)*3600.

    c = c + 1

# Calculate the IWP perturbations between each of the H[4:5] values.
IWP = np.loadtxt(basedir + '/output/IWP_iterator.txt')
dIWP = np.zeros(len(IWP)-1)
for j,iwp in enumerate(IWP):
    if j == 0:
       continue
    else:
       dIWP[j-1] = IWP[j] - IWP[j-1]

# Sensitivity of the LW/SW heating to Tt/Tb (1-K perturbations).
# Can be conceptualized moving from lower to upper atm or vice versa....
dHdT = np.zeros((6,H.shape[1]-1,H.shape[2]))
for j in np.arange(H.shape[1]-1,0,-1):
    print(j)
    dHdT[:4,j-1] = -1*(H[:4,j] - H[:4,j-1])
#for j in np.arange(H.shape[1]-1):
#    dHdT[:4,j] = (H[:4,j+1] - H[:4,j])
#    dHdT[4:,j] = (H[4:,j+1] - H[4:,j])/dIWP[j]

# x-tick labels become more complicated as some columns are zero.
xtlbl = [[] for i in np.arange(3)]

# The same columns are all 0 for the LW/SW fluxes from the Tt tests.
idx = np.argwhere(np.all(dHdT[0]==0,axis=1))
dHdT_Tb = np.zeros((2,dHdT.shape[1]-len(idx),dHdT.shape[2]))
dHdT_Tb[:2] = np.delete(dHdT[:2],idx,axis=1)
xtlbl[0] = tb[np.argwhere(~np.all(dHdT[0]==0,axis=1))[:,0]]

# The same columns are all 0 for the LW/SW fluxes from the Tb tests.
idx = np.argwhere(np.all(dHdT[2]==0,axis=1))
dHdT_Tt = np.zeros((2,dHdT.shape[1]-len(idx),dHdT.shape[2]))
dHdT_Tt[:2] = np.delete(dHdT[2:4],idx,axis=1)
xtlbl[1] = tt[np.argwhere(~np.all(dHdT[2]==0,axis=1))[:,0]]

# The third set of x-tick labels are IWP values. Factor of 1000 for g kg-1.
xtlbl[2] = [round(i*1000,2) for i in IWP]

for i in np.arange(6):
    print('min / max of sens matrix row ' + str(i) + ': ' + str(np.nanmin(dHdT[i])) + '  ' + str(np.nanmax(dHdT[i])))

fig, ax = plt.subplots(nrows=2,ncols=3,figsize=(11,11))

# Some formatting factors for the plot.. fontsize, tick intervals, and iterator.
fs = 13
c = 0
mm = np.array([[-1.1,1.1],[-1.1,1.1],[-2.5,7.5],[-0.1,0.1],[-0.1,0.1],[-5,7.5]])
lwsw = ['LW','SW']
xlbl = [r'$T_b$ with $T_t$ = 200 K',r'$T_t$ with $T_b$ = 237 K','IWP [g kg$^{-1}$']
lbl = [r'K day$^{-1}$ K$^{-1}$',r'K day$^{-1}$ K$^{-1}$',r'K day$^{-1}$ K$^{-1}$',\
       r'K day$^{-1}$ K$^{-1}$',r'K day$^{-1}$ per g kg$^{-1}$',r'K day$^{-1}$ per g kg$^{-1}$']

yt_tt_hl = np.array([int(t) for t in tt_hl[tropopause+3:49:5]])
field = [dHdT_Tb[0], dHdT_Tt[0], dHdT[4], dHdT_Tb[1], dHdT_Tt[1], dHdT[5]]

for i in np.arange(2):
    for j in np.arange(3):
        sns.heatmap(field[c].T,cmap=cm.bwr,ax=ax[i,j],norm=MidpointNormalize(midpoint=0.,vmin=mm[c,0],vmax=mm[c,1]),\
                   cbar_kws={'label':lbl[c]})

        # Adjust tick labels
        if j == 0:
           ax[i,j].set_ylabel('Temperature [K]',fontsize=fs)
        if i == 1:
           ax[i,j].set_xlabel(xlbl[j],fontsize=fs)

        ax[i,j].set_ylim([49,tropopause+3]) # ylim is done based on indices not values.
        ax[i,j].set_yticks(np.arange(tropopause+3,49,5))
        ax[i,j].set_yticklabels(yt_tt_hl,rotation=45,fontsize=fs-2)
        ax[i,j].set_xticks(np.arange(0,len(xtlbl[j]),2))  # Label every 2nd tick.
        ax[i,j].set_xticklabels(xtlbl[j][::2],rotation=45,fontsize=fs-2)

        # Where is 200 K?
        ax[i,j].plot([0,len(xtlbl[j])],[t200+tropopause,t200+tropopause],linewidth=1,color='k',linestyle='--')
        # Where is 237 K?
        ax[i,j].plot([0,len(xtlbl[j])],[t237+tropopause+3,t237+tropopause+3],linewidth=1,color='k',linestyle='--')

        c = c + 1

#fig.savefig('../figures/heating_matrix_qi2.pdf',bbox_inches='tight')
plt.show()
