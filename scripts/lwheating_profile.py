import matplotlib.pyplot as plt
import numpy as np
import os
from os.path import dirname

# Pressure at half levels
basedir = dirname(os.getcwd())

ellingson = np.genfromtxt(basedir + '/output/tropical_profile_ellingson_250m_formatted_top2bottom.txt')
pp_fl = ellingson[:,1]
pp_hl = np.zeros((82,))
for i in np.arange(1,81):
    pp_hl[i] = (pp_fl[i-1] + pp_fl[i])/2.0
pp_hl[0] = pp_fl[0] - (pp_hl[1] - pp_fl[0])
pp_hl[81] = pp_fl[80] + (pp_fl[80] - pp_hl[80])

# Heat capacity [J kg-1 K-1]
cp = 1.08*10**(3)
# Gravity [m s-2]
g = 9.8
# Temperatures at which qi was perturbed
#tt = np.asarray([261,262,264,266,267,269])
tt = np.asarray([202,204,205,207,209,210])

farbe = ['red','orange','gold','green','blue','purple']
#lbl = ['261 K','262 K','264 K','266 K','267 K','269 K']
lbl = ['202 K','204 K','205 K','207 K','209 K','210 K']

fs = 13
fig = plt.figure()
for c,tk in enumerate(tt):
    lw = np.genfromtxt(basedir + '/output/q1_3lev/lwflxatm-test' + str(tk) + '_q1.txt')
    # First column is all-sky. Second column is clear-sky. 
    # -1 is to define outgoing fluxes as positive.
    lwcld = -1.*(lw[0]-lw[1])

    # Calculate the longwave heating rate. 
    # Factor of 3600 converts from K s-1 to K day-1.
    H = g/cp*np.gradient(lwcld,pp_hl*100.)*3600
    plt.plot(H,pp_hl,color=farbe[c],linewidth=1.25,label=lbl[c])
    plt.ylabel('Pressure [hPa]',fontsize=fs)
    plt.xlabel(r'Cld LW heating rate [K day$^{-1}$]',fontsize=fs)

plt.legend()
plt.gca().invert_yaxis()
#fig.savefig('../figures/lwheating_qi1[2].pdf',bbox_inches='tight')
plt.show()
