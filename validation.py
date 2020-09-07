# This script takes a random timestep (hours after sim start) and gridcell, returns
# the array of inputs needed to run the offline RRTM, and prints the SW / LW fluxes
# generated by the online radiation codes.

def validation(whichhour, gridcell, testname):
    import numpy as np
    import xarray as xr
    import pandas as pd
    import sys

    # This was an attempt to run from command line more succintly.
    # gridcell = int(sys.argv[1])
    # whichhour = int(sys.argv[2])

    basedir0 = '/work/bb1131/b380873/tropic_vis/remapping/'
    #basedir1 = '/work/bb1131/b380873/tropic_run5_output/'
    basedir1 = '/scratch/b/b380873/tropic_run5/'
    basedir2 = '/work/bb1131/b380459/TROPIC/'

    # File <whichhour> of the 3D fields, file 2*<whichhour> of the 2D ones.
    n = whichhour
    m = 2*n

    while(len(str(n)) != 4):
         n = '0' + str(n)
    while(len(str(m)) != 4):
         m = '0' + str(m)

    # Read in the lat/lon associated with the gridcell.
    print('Read lat, lon, time, and SZA.')
    fi = xr.open_dataset(basedir0 + 'icon-gridcell_latlon_55e115e5s40n.nc')
    lat = fi.latitude.values[gridcell]
    #lon = fi.longitude.values[gridcell]
    del fi

    # Read in the time associated with timestep, lat/lon associated with gridcell.
    # Read in the cosine of the solar zenith angle and take its arccos.
    fi0 = xr.open_dataset(basedir1 + 'CLCONV_2D_icon_tropic_' + m + '.nc')
    time = pd.to_datetime(fi0.time.values[0])
    day = time.day
    mon = time.month
    hour = time.hour
    minute = time.minute
    cosmu0 = fi0.cosmu0.values[0,gridcell]
    sza = np.arccos(cosmu0)
    del fi0

    # There is a time dimension only for the aerosol fields.
    # Extract only the specificed gridcell.
    print('Read external parameters')
    fi1 = xr.open_dataset(basedir2 + '/extpar/extpar_icon-grid_tropic_55e115e5s40n_R2500m_bitmap.nc')
    zland = fi1.FR_LAND.values[gridcell]
    lwemiss = fi1.EMIS_RAD.values[gridcell]
    aer_bc = fi1.AER_BC.values[mon-1,gridcell]   # Zero-indexing in Python --> mon-1
    aer_du = fi1.AER_DUST.values[mon-1,gridcell]
    aer_org = fi1.AER_ORG.values[mon-1,gridcell]
    aer_so4 = fi1.AER_SO4.values[mon-1,gridcell]
    aer_ss = fi1.AER_SS.values[mon-1,gridcell]
    del fi1

    # Read in mass mixing ratios of vapor, liquid, and ice.
    print('Read 3D cloud fields.')
    fi2 = xr.open_dataset(basedir1 + 'CLCONV_3D_icon_tropic_' + n + '.nc')
    xm_vap = fi2.tot_qv_dia.values[0,:,gridcell]
    xm_ice = fi2.tot_qi_dia.values[0,:,gridcell]
    xm_liq = fi2.tot_qc_dia.values[0,:,gridcell]
    klev = len(fi2.tot_qc_dia.height)
    # Read in the CDNC and cloud fraction profiles.
    cdnc = fi2.acdnc.values[:,gridcell]
    cld_frc = fi2.clc.values[0,:,gridcell]
    del fi2

    # Read in the pressure and temperature at full and half levels.
    print('Read pressure and temperature.')
    fi3 = xr.open_dataset(basedir1 + 'WINDTH_3D_icon_tropic_' + n + '.nc')
    tk_fl = fi3.temp.values[0,:,gridcell]
    tk_hl = fi3.temp_ifc.values[0,:,gridcell]
    tk_sfc = fi3.t_g.values[0,gridcell]
    pp_fl = fi3.pres.values[0,:,gridcell]
    pp_hl = fi3.pres_ifc.values[0,:,gridcell]
    pp_sfc = fi3.pres_sfc.values[0,gridcell]
    del fi3

    # Read in surface albedos.
    print('Read surface albedos.')
    fi4 = xr.open_dataset(basedir1 + 'RAD_2D_icon_tropic_' + m + '.nc')
    alb_vis_dir = fi4.albvisdir.values[0,gridcell]
    alb_nir_dir = fi4.albnirdir.values[0,gridcell]
    alb_vis_dif = fi4.albvisdif.values[0,gridcell]
    alb_nir_dif = fi4.albnirdif.values[0,gridcell]
    del fi4

    # Stack all values together, flatten the list, and save in a txt file.
    inputs = [[klev], [day], [hour], [minute], [lat], [zland], [sza], [alb_vis_dir], [alb_nir_dir],\
              [alb_vis_dif], [alb_nir_dif], [lwemiss], list(pp_fl), list(pp_hl), [pp_sfc], list(tk_fl),\
              list(tk_hl), [tk_sfc], list(xm_vap), list(xm_liq), list(xm_ice), list(cdnc), list(cld_frc),\
              [aer_ss], [aer_org], [aer_bc], [aer_so4]]

    # The name of the text file is based on an input string.
    with open('rrtm_inputs_' + testname + '.txt','w') as f:
         for val in inputs:
             for i,elem in enumerate(val):
                if i != len(val)-1:
                   f.write(str(elem) + ' ')
                else:
                   f.write(str(elem))
             f.write('\n')
    f.close()

    # Print the flux values against which we check.
    fi5 = xr.open_dataset(basedir1 + 'RAD_3D_icon_tropic_' + n + '.nc')
    swflxall = fi5.swflxall.values[0,:,gridcell]
    swflxclr = fi5.swflxclr.values[0,:,gridcell]
    lwflxall = fi5.lwflxall.values[0,:,gridcell]
    lwflxclr = fi5.lwflxclr.values[0,:,gridcell]

#if __name__ == '__main__':
#  validation()
