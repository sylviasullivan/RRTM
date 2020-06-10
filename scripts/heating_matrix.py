import matplotlib.pyplot as plt
import matplotlib.colors as colors
import os
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
# Temperatures at which qi was perturbed
#tt = np.asarray([196,197,199,200,202,204,205,207,209,210,212,214,215,217,\
#        218,220,222,225,227,229,230,232,234,235,237,239,240,242,244,\
#         245,247,249,251,252,254,256,257,260,261,262,264,266,267,269])
tt = np.asarray([202,204,205,207,209,210,212,214,215,217,\
        218,220,222,225,227,229,230,232,234,235,237,239,240,242,244,\
         245,247,249,251,252,254,256,257,260,261])

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

fs = 13

# Matrix to hold the radiative heating profiles
H = np.zeros((4,tt.shape[0],lev+1))

for c,tk in enumerate(tt):
    lw1 = np.genfromtxt(basedir + '/output/q2_1lev/lwflxatm-test' + str(tk) + '_q2.txt')
    sw1 = np.genfromtxt(basedir + '/output/q2_1lev/swflxatm-test' + str(tk) + '_q2.txt')
    lw5 = np.genfromtxt(basedir + '/output/q2_3lev/lwflxatm-test' + str(tk) + '_q2.txt')
    sw5 = np.genfromtxt(basedir + '/output/q2_3lev/swflxatm-test' + str(tk) + '_q2.txt')

    # First column is all-sky. Second column is clear-sky.
    # The factor -1 is to define outgoing fluxes as positive.
    lwcld1 = -1.*(lw1[0]-lw1[1])
    swcld1 = -1.*(sw1[0]-sw1[1])
    lwcld5 = -1.*(lw5[0]-lw5[1])
    swcld5 = -1.*(sw5[0]-sw5[1])

    # Calculate the longwave and shortwave heating rates.
    # Factor of 3600 converts from K s-1 to K day-1.
    H[0,c] = g/cp*np.gradient(lwcld1,pp_hl*100.)*3600.
    H[1,c] = g/cp*np.gradient(swcld1,pp_hl*100.)*3600.
    H[2,c] = g/cp*np.gradient(lwcld5,pp_hl*100.)*3600.
    H[3,c] = g/cp*np.gradient(swcld5,pp_hl*100.)*3600.


fig, ax = plt.subplots(nrows=2,ncols=2)
# Some factors to make the plot.. fontsize, tick intervals, and iterator.
fs = 13
interval = 5
c = 0
mm = [[-1.1,1.7],[-0.11,0.17],[-1.1,1.7],[-0.11,0.17]]

yt_tt_hl = np.array([int(t) for t in tt_hl[::interval]])
yt_pp_hl = np.array([int(t) for t in pp_hl[::interval]])
xt_tt = np.array([int(t) for t in tt[::interval]])

for i in np.arange(2):
    for j in np.arange(2):
        sns.heatmap((H[c].T),norm=colors.SymLogNorm(vmin=H[c].min(),vmax=H[c].max(),linthresh=0.001),
            cmap=cm.bwr,xticklabels=xt_tt,yticklabels=yt_tt_hl,cbar_kws={'label':r'K day$^{-1}$'},ax=ax[i,j])

        # Where is the melting layer?
        ax[i,j].plot([0,tt.shape[0]],[melting_layer,melting_layer],linewidth=1,color='k',linestyle='--')
        # Where is the tropopause?
        ax[i,j].plot([0,tt.shape[0]],[tropopause,tropopause],linewidth=1,color='k',linestyle='--')
        # Adjust tick labels
        ax[i,j].set_xticks(np.arange(0,tt.shape[0],interval))
        ax[i,j].set_xticklabels(ax[i,j].get_xticklabels(),rotation=45,fontsize=fs-2)
        ax[i,j].set_yticks(np.arange(0,tt_hl.shape[0],interval))
        ax[i,j].set_yticklabels(ax[i,j].get_yticklabels(),rotation=45,fontsize=fs-2)
        #ax[i,j].set_ylabel('Temperature level of LW heating perturbation [K]',fontsize=fs)
        #ax[i,j].set_xlabel('Temperature level of $q_i$ perturbation [K]',fontsize=fs)
        c += 1


#plt.ylabel('Temperature level of SW heating perturbation [K]',fontsize=fs)
#plt.xlabel('Temperature level of $q_i$ perturbation [K]',fontsize=fs)

#fig.savefig('../figures/lwheating_qi002_matrix.pdf',bbox_inches='tight')
plt.show()
