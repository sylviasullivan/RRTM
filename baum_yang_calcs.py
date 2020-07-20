import sys
import numpy as np
import xarray as xr
from netCDF4 import Dataset

# Below values are taken from W. Zhao Atm. Res. 204 (2018) 37-53
# Spectral ranges:
# SW - (3.077,3.846) (2.5-3.077) (2.15-2.5) (1.942-2.15) (1.626-1.942) (1.299-1.626) (1.242-1.299)
# (0.778-1.242) (0.625-0.778) (0.442-0.625) (0.345-0.442) (0.263-0.345) (0.2-0.263) (3.846-12.195)



# SW mass extinction coefficient
# beta = IWC (a0 + a1 / De + a2 / De^2)
# a0 [m-1 / (g m-3)], a1 [1 / (g m-3)], a2 [m / (g m-3)]
a = np.array([[-0.301e-2,-0.2379e-2,-0.1787e-2,-0.1864e-2,-0.1819e-2,-0.1640e-2,-0.1836e-2,\
              -0.1695e-2,-0.1574e-2,-0.135e-2,-0.1188e-2,-0.1088e-2,-0.9619e-3,-0.8349e-2],\
             [0.3318e1,0.3208e1,0.3170e1,0.3183e1,0.3186e1,0.3171e1,0.3193e1,0.3173e1,0.3154e1,\
              0.3131e1,0.3115e1,0.3106e1,0.3087e1,0.3829e1],\
             [0.1084e1,-0.4263e1,0.1997e1,0.1596e1,0.6773,0.5581,0.4244,0.2865,\
              0.1543,0.7913e-3,-0.2874,-0.4073,-0.6622,-0.683e1]])

# SW coalbedo
# 1 - omega = b0 + b1 * De + b2 * De^2 + b3 * De^3
b = np.array([[0.2693,0.2075,-0.2841e-3,-0.2892e-2,-0.1347e-4,-0.1188e-4,0.359e-4,0.4152e-4,\
               0,0,0,0,0,1.715],\
              [0.4137e-2,0.8103e-3,0.9704e-3,0.3478e-2,0.7488e-3,0.8368e-3,0.6289e-4,0.7141e-5,\
               0,0,0,0,0,0.695e-3],\
              [-0.3817e-4,0.5698e-5,-0.8236e-7,-0.1556e-4,0.1627e-6,-0.9913e-6,0.1237e-6,\
               0.3644e-7,0,0,0,0,0,-0.6358e-4],\
              [0.1230e-6,-0.4896e-7,-0.1237e-7,0.1665e-7,-0.9766e-8,-0.7163e-8,-0.738e-9,\
               -0.1417e-9,0,0,0,0,0,0.1991e-6]])

# SW asymmetry factor
# g = c0 + c1 * De + c2 * De^2 + c3 * De^3
c = np.array([[0.8011,0.8837,0.8172,0.7931,0.7811,0.7745,0.7717,0.7681,0.7662,0.7646,0.7615,\
               0.7562,0.7358,0.8222],\
              [0.4177e-2,0.1185e-2,0.1637e-3,0.1563e-2,0.7757e-3,0.9057e-3,0.6579e-3,0.7369e-3,\
               0.7545e-3,0.7328e-3,0.6768e-3,0.6168e-3,0.4968e-3,0.3681e-2],\
              [-0.4192e-4,-0.5850e-5,0.1079e-4,-0.1117e-5,0.1888e-5,-0.6592e-7,-0.6515e-6,\
               -0.2626e-5,-0.3408e-5,-0.3405e-5,-0.2815e-5,-0.1936e-5,0.1301e-6,-0.3126e-4],\
              [0.1422e-6,0.3932e-8,-0.6696e-7,-0.331e-7,-0.2739e-7,-0.1941e-7,-0.1075e-7,\
               -0.7498e-9,0.3847e-8,0.4707e-8,0.285e-8,-0.8715e-9,-0.1025e-7,0.8696e-7]])

# LW absorption coefficient
# beta_a = IWC (d0 + d1 / De + d2 / De^2 + d3 / De^3)
# d0 [m-1 / (g m-3)], d1 [1 / (g m-3)], d2 [m / (g m-3)], d3 [m2 / (g m-3)]
d = np.array([[-0.2214e-2,0.134e-2,-0.2015e-2,-0.1185e-2,-0.1004e-2,-0.6071e-3,0.8331e-3,\
               0.9738e-3,0.5614e-3,0.4112e-4,-0.1463e-4,0.1801e-2,0.7473e-3,0.5736e-3,0.1344e-2,\
               0.1228e-2],\
              [0.1813e1,0.1030e1,0.1675e1,0.1513e1,0.1455e1,0.1444e1,0.1354e1,0.1318e1,0.1374e1,\
               0.1445e1,0.1449e1,0.1149e1,0.1343e1,0.1381e1,0.1239e1,0.1192e1],\
              [-0.1587e2,-0.1509e2,-0.1032e2,0.4817e1,0.7128e1,0.8203,-0.1143e2,-0.1087e2,\
               -0.8818e1,-0.5474e1,-0.4716e1,-0.1191e2,-0.9834e1,-0.9015e1,-0.1174e2,-0.8474e1],\
              [0.492e2,0.699e2,0.1224e2,-0.6885e2,-0.6147e2,-0.2737e2,0.3563e2,0.3283e2,0.1833e2,\
               -0.2172e1,-0.4387e1,0.4768e2,0.3044e2,0.2439e2,0.411e2,0.3314e2]])

# How many shortwave / longwave bands?
SW_band = 14
LW_band = 16

# Range of equivalent radii [m]
cldopt = xr.open_dataset('rrtm_cldopt.nc')
re = cldopt['re_crystal'].values*10**(-6)

# Range of wavelengths [m]
lamb = cldopt['wavelength'].values*10**(-6)

# Matrices of radii and longwave / shortwave wavelength values.
DeLW_M, _ = np.meshgrid(re*2,lamb[:LW_band])
DeSW_M, _ = np.meshgrid(re*2,lamb[LW_band:])

# SW mass extinction coefficient divided by IWC
# (Multiply by IWC elsewhere to return mass extinction coefficient.)
beta_per_IWC = np.array([a[0],]*len(re)).T + np.array([a[1],]*61).T/DeSW_M + \
               np.array([a[2],]*len(re)).T/DeSW_M**2

# SW co-albedo
omega = 1 - (np.array([b[0],]*len(re)).T + np.array([b[1],]*len(re)).T*DeSW_M + \
        np.array([b[2],]*len(re)).T*DeSW_M**2 + np.array([b[3],]*len(re)).T*DeSW_M**3)

# SW asymmetry factor
g = np.array([c[0],]*len(re)).T + np.array([c[1],]*len(re)).T*DeSW_M + \
    np.array([c[2],]*len(re)).T*DeSW_M**2 + np.array([c[3],]*61).T*DeSW_M**3
f = 1/(2*omega)

# LW mass absorption coefficient dividied by IWC
# (Multiply by IWC elsewhere to return mass absorption coefficient.)
betaa_per_IWC = np.array([d[0],]*len(re)).T + np.array([d[1],]*len(re)).T/DeLW_M + \
                np.array([d[2],]*len(re)).T/DeLW_M**2 + np.array([d[3],]*len(re)).T/DeLW_M**3

# Generate a netcdf file with all the optical properties according to Baum-Yang.
# Create dimensions of equivalent radius and wavelength.
by_nc = Dataset('baum_yang_cldopt.nc','w',format='NETCDF4')
req = by_nc.createDimension('re_crystal',re.shape[0])
wavel = by_nc.createDimension('wavelength',lamb.shape[0])

# Some file descriptors
by_nc.description = "Cloud optical properties for ice clouds ONLY at 30 wavelengths"
by_nc.author = "Sylvia Sullivan (sylvia.sullivan@kit.edu) following Zhao et al. AR (204) 37-53, 2018"
by_nc.title = "Baum-Yang Cloud Optical Properties"

# Create variables with appropriate dimensions.
# These are formatted identically to those in the existing rrtm_cldopt.nc
wavel_vals  = by_nc.createVariable('wavelength',np.float32,('wavelength'))
waven_vals  = by_nc.createVariable('wavenumber',np.float32,('wavelength'))
re_crystal  = by_nc.createVariable('re_crystal',np.float32,('re_crystal'))
rrtm_band   = by_nc.createVariable('rrtm_band',np.float32,('wavelength'))

eps_crystal = by_nc.createVariable('extinction_per_mass_crystal',np.float32,\
              ('wavelength','re_crystal'))
alb_crystal = by_nc.createVariable('co_albedo_crystal',np.float32,\
              ('wavelength','re_crystal'))
asy_crystal = by_nc.createVariable('asymmetry_factor_crystal',np.float32,\
              ('wavelength','re_crystal'))

# Set the values of the variables.
wavel_vals[:] = lamb * 10**(6) # Convert back to um.
wavel_vals.long_name = 'Wavelength (at center of band)'
wavel_vals.units = 'microns'
wavel_vals.fill_value = -9999.

waven_vals[:] = 1 / lamb * 10**(-2)
waven_vals.long_name = 'Wavenumber (at mid wavelength)'
waven_vals.units = 'cm-1'
waven_vals.fill_value = -9999.

re_crystal[:] = re * 10**6  # Convert back to um.
re_crystal.long_name = 'Effective Radius (Crystals)'
re_crystal.units = 'microns'
re_crystal.fill_value = -9999.

#rrtm_band[:] =
rrtm_band.long_name = 'Band number in RRTMG Code'
rrtm_band.fill_value = -999.

eps_crystal[:] = np.concatenate((beta_per_IWC,betaa_per_IWC),axis=0)
eps_crystal.long_name = 'Extinction Coefficient [non-spherical ice]'
eps_crystal.units = 'm-1'
eps_crystal.fill_value = -9999.

# Zeros have to be appended in the LW for the SW-only properties.
omega_2 = np.concatenate((np.zeros((LW_band,len(re))),1-omega),axis=0)
alb_crystal[:] = omega_2
alb_crystal.long_name = 'Co-Albedo (1-ssa), [non-spherical ice]'
alb_crystal.units = '1'
alb_crystal.fill_value = -9999.

g_2 = np.concatenate((np.zeros((LW_band,len(re))),g),axis=0)
asy_crystal[:] = g_2
asy_crystal.long_name = 'Asymmetry factor, [non-spherical ice]'
asy_crystal.units = '1'
asy_crystal.fill_value = -9999.

