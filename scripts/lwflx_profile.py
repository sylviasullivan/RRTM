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

print(pp_hl)


farbe = ['red','orange','gold','green','blue','purple']
tt = np.asarray([261,262,264,266,267,269])

fig = plt.figure()
for c,tk in enumerate(tt):
    lw = np.genfromtxt(basedir + '/output/q1_1lev/lwflxatm-test' + str(tk) + '_q1.txt')
    # First column is all-sky. Second column is clear-sky. 
    lwcld = lw[0]-lw[1]
    plt.plot(lwcld,pp_hl,color=farbe[c],linewidth=1.25)

plt.gca().invert_yaxis()
plt.show()
