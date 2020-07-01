import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.colors as colors

h2o = xr.open_dataset('/home/sylvia/Documents/rrtm/h2o_ctm_Ts300_rh0.75_gamma7.nc')
basedir = '/home/sylvia/Documents/rrtm/figures/'
# Gravitational acceleration [m s-2]
g = 9.8
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
# Ice density [kg m-3]
rho_ice = 920


ice_opt = np.genfromtxt('/home/sylvia/Documents/rrtm/output/ice_optical_properties.txt')
# Wavelength [um]
lamb = ice_opt[:,0]
# Wavenumbers [cm-1] as above
wn_ice = 1./lamb*10**(4)
# Imaginary part of the index of refraction
mim = ice_opt[:,2]
# Calculate the ice mass absorption coefficient from m_im. [m2 kg-1]
# Note that we have to make an assumption about ice density here.
sigma_ice = 4*np.pi*mim/(lamb*10**(-6))/rho_ice

# Plot the absorption coefficients versus wavenumber for a temperature = 220 K.
fig2, ax = plt.subplots(nrows=1,ncols=2,figsize=(11,5))
fs = 15
sigma_wv_T220 = sigma_wv[np.argmin(T-220)]
print('Min abs coeff at 220 K: ' + str(np.nanmin(sigma_wv_T220)) + ' m2 kg-1, Max abs coeff at 220 K: ' + 
        str(np.nanmax(sigma_wv_T220)) + ' m2 kg-1')
ax[0].plot(wn_wv,sigma_wv_T220)
ax[0].set_xlabel(r'Wavenumber [cm$^{-1}$]',fontsize=fs)
ax[0].set_ylabel(r'$\kappa_{wv}$ [m$^{2}$ kg$^{-1}$]',fontsize=fs)
ax[0].set_yscale('log')
ax[0].tick_params(labelsize=fs,axis='both',which='major')
ax[0].set_xlim([150,1450])
ax[0].set_ylim([10**(-8),100])
ax[0].tick_params(labelsize=fs)

ax[1].plot(wn_ice,sigma_ice)
ax[1].set_xlabel(r'Wavenumber [cm$^{-1}$]',fontsize=fs)
ax[1].set_ylabel(r'$\kappa_{ice}$ [m$^{2}$ kg$^{-1}$]',fontsize=fs)
ax[1].set_yscale('log')
ax[1].tick_params(labelsize=fs,axis='both',which='major')
ax[1].set_xlim([150,1450])
ax[1].set_ylim([1,1000])
ax[1].tick_params(labelsize=fs)

#fig2.savefig(basedir + 'wv_ir_absorption.pdf')
#plt.show()

# Plot the emission spectrum of ice
tau = -1.*np.integral(sigma_ice*2*10**(-5)/g

fig3 = plt.figure()
plt.plot(lamb,1-np.e**(-1.*tau)
plt.ylabel(r'$\kappa_{ice}$ [m$^{2}$ kg$^{-1}$]',fontsize=fs)
plt.gca().set_yscale('log')
#plt.xlabel(r'Wavenumber [cm$^{-1}$]',fontsize=fs)
plt.xlabel(r'Wavelength [$\mu$m]',fontsize=fs)
plt.xlim([0.1,1000])
plt.gca().set_xscale('log')


