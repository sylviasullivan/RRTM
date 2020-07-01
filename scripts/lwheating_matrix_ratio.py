import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import seaborn as sns
from matplotlib import cm
import os
from os.path import dirname

exec(open('MidpointNormalize.py').read())
basedir = dirname(os.getcwd())

# Number of levels over which the RRTM calculations are done
lev = 81
# Heat capacity [J kg-1 K-1]
cp = 1.08*10**(3)
# Gravity [m s-2]
g = 9.8
# Temperatures at which qi was perturbed
tt = np.asarray([196,197,199,200,202,204,205,207,209,210,212,214,215,217,\
        218,220,222,225,227,229,230,232,234,235,237,239,240,242,244,\
         245,247,249,251,252,254,256,257,260,261,262,264,266,267,269])

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
# Matrix to hold the LW radiative heating profiles
Hmat1 = np.zeros((tt.shape[0],lev+1))
Hmat2 = np.zeros((tt.shape[0],lev+1))

for c,tk in enumerate(tt):
    lw1 = np.genfromtxt(basedir + '/output/q1_1lev/lwflxatm-test' + str(tk) + '_q1.txt')
    # First column is all-sky. Second column is clear-sky.
    lwcld1 = lw1[0]-lw1[1]
    
    lw2 = np.genfromtxt(basedir + '/output/q2_1lev/lwflxatm-test' + str(tk) + '_q2.txt')
    lwcld2 = lw2[0] - lw2[1]

    # Calculate the longwave heating rate.
    # Factor of 3600 converts from K s-1 to K day-1.
    H1 = g/cp*np.gradient(lwcld1,pp_hl*100.)*3600
    H2 = g/cp*np.gradient(lwcld2,pp_hl*100.)*3600

    Hmat1[c] = H1
    Hmat2[c] = H2

Hratio = (Hmat2 - Hmat1) / Hmat1
print(np.nanmin(Hratio),np.nanmean(Hratio),np.nanmax(Hratio))

fig = plt.figure()
fs = 13
interval = 5
yt_tt_hl = np.array([int(t) for t in tt_hl[::interval]])
yt_pp_hl = np.array([int(t) for t in pp_hl[::interval]])
xt_tt = np.array([int(t) for t in tt[::interval]])

#ax = sns.heatmap(np.transpose(Hmat),norm=MidpointNormalize(midpoint=0,vmin=-0.25,vmax=0.3),
ax = sns.heatmap(np.transpose(Hratio),norm=MidpointNormalize(midpoint=0,vmin=-10,vmax=100),#norm=colors.SymLogNorm(vmin=-10,vmax=100,linthresh=0.1),
        #norm=colors.SymLogNorm(Hratio.min(),vmax=Hratio.max(),linthresh=0.01),
        cmap=cm.bwr,xticklabels=xt_tt,yticklabels=yt_tt_hl,cbar_kws={'label':r'$q_{i,2}$ / $q_{i,1}$'})

# Where is the melting layer?
ax.plot([0,tt.shape[0]],[melting_layer,melting_layer],linewidth=1,color='k',linestyle='--')
# Where is the tropopause?
ax.plot([0,tt.shape[0]],[tropopause,tropopause],linewidth=1,color='k',linestyle='--')
# Adjust tick labels
ax.set_xticks(np.arange(0,tt.shape[0],interval))
ax.set_xticklabels(ax.get_xticklabels(),rotation=45,fontsize=fs-2)
ax.set_yticks(np.arange(0,tt_hl.shape[0],interval))
ax.set_yticklabels(ax.get_yticklabels(),rotation=45,fontsize=fs-2)

plt.ylabel('Temperature level of LW heating perturbation [K]',fontsize=fs)
plt.xlabel('Temperature level of $q_i$ perturbation [K]',fontsize=fs)
#fig.savefig('../figures/lwheating_qi1_matrix.pdf',bbox_inches='tight')
plt.show()
