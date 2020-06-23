#!/bin/bash
module load gcc/7.1.0; export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/sw/rhel6-x64/netcdf/netcdf_fortran-4.4.3-gcc71/lib

#gfortran -g -Wall -Wmaybe-uninitialized -ffpe-trap=zero,overflow,underflow -fbounds-check -fbacktrace -fmax-errors=10 -I/sw/rhel6-x64/netcdf/netcdf_fortran-4.4.3-gcc71/include -L/sw/rhel6-x64/netcdf/netcdf_fortran-4.4.3-gcc71/lib -lnetcdff mo_lrtm_kgs.f90 mo_lrtm_setup.f90 mo_lrtm_rtrnmr.f90 mo_lrtm_netcdf.f90 mo_lrtm_coeffs.f90 mo_lrtm_taumol.f90 o3_gems_data.f90 o3_gems.f90 mo_lrtm.f90 mo_srtm_kgs.f90 srtm_kgb16.f90 srtm_kgb17.f90 srtm_kgb18.f90 srtm_kgb19.f90 srtm_kgb20.f90 srtm_kgb21.f90 srtm_kgb22.f90 srtm_kgb23.f90 srtm_kgb24.f90 srtm_kgb25.f90 srtm_kgb26.f90 srtm_kgb27.f90 srtm_kgb28.f90 srtm_kgb29.f90 mo_srtm_config.f90 mo_srtm_taumol.f90 mo_srtm.f90 rrtm.f90 rrtm_driver.f90

#LIBRARY AND INCLUDE - LOCAL
#-I/home/sylvia/libs/include -I/usr/include
#-L/home/sylvia/libs/lib -L/usr/lib/x86_64-linux-gnu

#LIBRARY AND INCLUE - DKRZ
#-I/sw/rhel6-x64/netcdf/netcdf_fortran-4.4.3-gcc71/include
#-L/sw/rhel6-x64/netcdf/netcdf_fortran-4.4.3-gcc71/lib
