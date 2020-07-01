import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import os
from os.path import dirname

basedir = dirname(os.getcwd())
h2o = xr.open_dataset(basedir + '/h2o_ctm_Ts300_rh0.75_gamma7.nc')

# Pressure levels [Pa]
p = h2o.p
# Wavenumbers [m-1], Factor of 100 to convert to [cm-1]
wn_wv = h2o.k/100.
# Temperature levels [K]
T = h2o.tabs
# Avogadros number [molecule mol-1]
Nav = 6.022*10**(23)
# Molar mass [kg mol-1]
Mh2o = 18.01528*10**(-3)
# Absorption coefficients [m2 molecule-1]
sigma_wv = h2o.sigma_h2o
sigma_wv = sigma_wv*Nav*Mh2o   # [m2 kg-1]
# Longwave cooling [K day-1 cm-1]
H = h2o.coo

print('Min pressure: ' + str(np.nanmin(p)) + ' Pa, Max pressure: ' + str(np.nanmax(p)) + 
      ' Pa, Pressure dims: ' + str(p.shape))
print('Min wv wavenumber: ' + str(np.nanmin(wn_wv)) + ' cm-1, Max wv wavenumber: ' + str(np.nanmax(wn_wv)) + 
      ' cm-1, wv wavenumber dims: ' + str(wn_wv.shape))
print('Min LW H: ' + str(np.nanmin(H)) + ' K day-1 cm-1, Max LW H: ' + str(np.nanmax(H)) + 
      ' K day-1 cm-1, wv abs coeff dims: ' + str(H.shape))

# Plot the space of absorption coefficients versus wavenumber and pressure / temperature.
fig = plt.figure()
plt.contourf(wn_wv,T,H,levels=20,vmin=-0.015,vmax=0.015) #,norm=colors.LogNorm(vmin=sigma_wv.min(),vmax=sigma_wv.max()))
plt.colorbar()
#fig.savefig(basedir + 'lwheating_T_wavenumber.pdf',bbox_inches='tight')

plt.show()
