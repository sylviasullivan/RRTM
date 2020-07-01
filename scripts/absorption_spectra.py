import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.colors as colors

h2o = xr.open_dataset('/home/sylvia/Documents/rrtm/h2o_ctm_Ts300_rh0.75_gamma7.nc')
basedir = '/home/sylvia/Documents/rrtm/figures/'
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

print('Min pressure: ' + str(np.nanmin(p)) + ' Pa, Max pressure: ' + str(np.nanmax(p)) + 
      ' Pa, Pressure dims: ' + str(p.shape))
print('Min wv wavenumber: ' + str(np.nanmin(wn_wv)) + ' cm-1, Max wv wavenumber: ' + str(np.nanmax(wn_wv)) + 
      ' cm-1, wv wavenumber dims: ' + str(wn_wv.shape))
print('Min wv abs coeff: ' + str(np.nanmin(sigma_wv)) + ' m2 kg-1, Max wv abs coeff: ' + str(np.nanmax(sigma_wv)) + 
      ' m2 kg-1, wv abs coeff dims: ' + str(sigma_wv.shape))

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

print('Min ice wavenumber: ' + str(np.nanmin(wn_ice)) + ' cm-1, Max ice wavenumber: ' + str(np.nanmax(wn_ice)) +
      ' cm-1, Wavenumber dims: ' + str(wn_ice.shape))
print('Min ice abs coeff: ' + str(np.nanmin(sigma_ice)) + ' m2 kg-1, Max wv abs coeff: ' + str(np.nanmax(sigma_ice)) +
      ' m2 kg-1, wv abs coeff dims: ' + str(sigma_ice.shape))

# Plot the space of absorption coefficients versus wavenumber and pressure / temperature.
fig = plt.figure()
plt.contourf(wn_wv,T,sigma_wv,levels=20,norm=colors.LogNorm(vmin=sigma_wv.min(),vmax=sigma_wv.max()))
plt.colorbar()
#fig.savefig(basedir + 'kwv_T_wavenumber.pdf',bbox_inches='tight')

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

# Plot the full range of ice absorption coefficients
fig3 = plt.figure()
plt.plot(lamb,sigma_ice)
plt.ylabel(r'$\kappa_{ice}$ [m$^{2}$ kg$^{-1}$]',fontsize=fs)
plt.gca().set_yscale('log')
#plt.xlabel(r'Wavenumber [cm$^{-1}$]',fontsize=fs)
plt.xlabel(r'Wavelength [$\mu$m]',fontsize=fs)
plt.xlim([0.1,1000])
plt.gca().set_xscale('log')

# Plot the Imaginary index of refraction versus wavelength as in Warren and Brandt 2008 Fig 4
fig4 = plt.figure()
plt.plot(lamb,mim)
plt.xlabel(r'Wavenumber [cm$^{-1}$]',fontsize=fs)
plt.xlim([5,20])
plt.ylim([10**(-2),1])
plt.gca().set_yscale('log')
plt.gca().tick_params(labelsize=fs,axis='both',which='major')
plt.show()

