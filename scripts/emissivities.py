import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib import colors
import sys

epsice = np.linspace(0.001,0.999,50)
epswv  = np.linspace(0.001,0.999,50)
basedir = '/home/sylvia/Documents/rrtm/figures/'
exec(open('MidpointNormalize.py').read())

# Create a matrix of results.
xx,yy = np.meshgrid(epsice,epswv)
zz = 0.5*xx/(1-xx) - yy

fs = 13
fig = plt.figure()
plt.plot(epswv,0.5*epsice/(1-epsice),linewidth=1.25,color='red',label=r'0.5$\varepsilon_{ice}$/(1-$\varepsilon_{ice}$)')
plt.plot(epswv,epswv,linewidth=1.25,color='blue',label=r'$\varepsilon_{wv}$')
plt.plot([0.5,0.5],[0,0.5],linewidth=0.75,linestyle='--',color='black')
plt.xlim([0,1])
plt.xlabel(r'$\varepsilon$',fontsize=fs)
plt.ylabel('Expression values',fontsize=fs-1)
plt.legend(fontsize=fs)
plt.gca().set_yscale('log')
plt.gca().tick_params(labelsize=fs)
fig.savefig(basedir + 'epsilons.pdf',bbox_inches='tight')

fig2 = plt.figure(figsize=(6,5))
ax = fig2.gca(projection='3d')
# Impose a log scale in z since this option does not exist for Axes3D.
#zz[zz > 0] = np.log10(zz[zz > 0])
#zz[zz < 0] = -1.*np.log10(np.abs(zz[zz < 0]))

Elim = 2
zz[zz > Elim] = np.nan
ax.plot_surface(xx,yy,np.zeros(xx.shape),color='gray',alpha=0.5)
ax.plot_surface(xx,yy,zz,cmap=cm.coolwarm,norm=MidpointNormalize(midpoint=0,vmin=-1,vmax=Elim),
        label=r'0.5$\varepsilon_{ice}$/(1-$\varepsilon_{ice}$)-$\varepsilon_{wv}$')
ax.set_xlabel(r'$\varepsilon_{ice}$',fontsize=fs)
ax.set_ylabel(r'$\varepsilon_{wv}$',fontsize=fs)
ax.set_zlabel(r'0.5$\varepsilon_{ice}$/(1-$\varepsilon_{ice}$)-$\varepsilon_{wv}$',fontsize=fs)
ax.xaxis.set_tick_params(labelsize=fs+1)
ax.yaxis.set_tick_params(labelsize=fs+1)
ax.zaxis.set_tick_params(labelsize=fs+1)
ax.view_init(elev=20,azim=48)
#fig2.savefig(basedir + 'epsilon_surf.pdf',bbox_inches='tight')
plt.show()
