# RRTM - Rapid Radiative Transfer Module

Offline radiative calculations. RRTM uses k-distributions from line-by-line calculations. 16 longwave bands (interval 1: 10-250 cm-1, interval 2: 250-500 cm-1, interval 3: 500-630 cm-1, interval 4: 630-700 cm-1 .. up to 3000 cm-1) and 14 shortwave bands (interval 1: 2600-3250 cm-1, interval 2: 3250-4000 cm-1, interval 3: 4000-4650 cm-1, interval 4: 4650-5150 cm-1 ... up to 50000 cm-1) with 140 quadrature points are used.

**rrtm**
Here we set trace gase species other than O3 (water vapor, nitrous oxide, CH4, halocarbons) as well as aerosol species (BC, organics, sea salt, dust).

**mo_lrtm_rtrnmr**
Upwelling and downwelling LW fluxes and heating rates are calculated from the input Planck function, thermodynamic profiles, and cloud fractions.

**mo_lrtm_taumol**
16 subroutines are used to calculate optical depth and Planck fractions per g-value and layer across the 16 LW bands.

**mo_lrtm_coeff**
Temperature and pressure interpolations for optical depths.

**mo_lrtm_netcdf** + **mo_lrtm_setup**
Read in and setup cloud optical properties, i.e. liquid droplet and ice crystal mass extinction coefficients, albedos, and asymmetry factors as a function of wavelength and effective radius

**mo_srtm**
The solar absorption spectrum and surface albedo are read in here.
*srtm_setcoef* - temperature and pressure dependence of SW optical depths
*srtm_spcvrt* - 2-stream SW flux calculations
*srtm_vrtqdr* - vertical quadrature
*srtm_reftra* - reflectivity and transmissivity of a given layer

**mo_srtm_taumol**
14 subroutines are used to calculate optical depth and Planck fractions per g-value and layer across the 14 SW bands.

**mo_srtm_kgs**
Absolute absorption coefficients across the 14 bands.

**o3_gems**
Geostationary Environment Monitoring Spectrometer (GEMS) climatological ozone fields

