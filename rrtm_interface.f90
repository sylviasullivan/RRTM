  SUBROUTINE rrtm_interface(                                              &
    ! input
    & klev            ,zland           ,pmu0            ,alb_vis_dir     ,&
    & alb_nir_dir     ,alb_vis_dif     ,alb_nir_dif     ,emis_rad        ,&
    & pp_fl           ,pp_hl           ,pp_sfc                           ,&
    & tk_fl           ,tk_hl           ,tk_sfc                           ,&
    & xm_vap          ,xm_liq          ,xm_ice                           ,&
    & cdnc            ,cld_frc         ,xm_o3           ,aotss           ,&
    & aotorg          ,aotbc           ,aotso4          ,aotdu           ,&
    ! >> sylvia_20200305
    ! Changing inputs zaeq1,zaeq2,zaeq3,zaeq4,zaeq5 to aot* variables.
    ! output
    & flx_lw_net      ,flx_sw_net      ,flx_lw_net_clr  ,flx_sw_net_clr  ,&
    & flx_uplw_sfc    ,flx_upsw_sfc    ,flx_uplw_sfc_clr,flx_upsw_sfc_clr,&
    ! optional output
    & flx_upsw_toa                                                        )
    
    INTEGER,INTENT(in)  ::         &
      &  klev                               !< number of levels

    REAL(wp),INTENT(in) ::         &
      &  zland,                    & !< land-sea mask. (1. = land, 0. = sea/lakes)
      &  pmu0,                     & !< mu0 for solar zenith angle
      &  alb_vis_dir,              & !< surface albedo for vis range and dir light
      &  alb_nir_dir,              & !< surface albedo for NIR range and dir light
      &  alb_vis_dif,              & !< surface albedo for vis range and dif light
      &  alb_nir_dif,              & !< surface albedo for NIR range and dif light
      &  emis_rad,                 & !< longwave surface emissivity
      &  pp_fl(klev),              & !< full level pressure in Pa
      &  pp_hl(klev+1),            & !< half level pressure in Pa
      &  pp_sfc,                   & !< surface pressure in Pa
      &  tk_fl(klev),              & !< full level temperature in K
      &  tk_hl(klev+1),            & !< half level temperature in K
      &  tk_sfc,                   & !< surface temperature in K
      &  xm_vap(klev),             & !< specific humidity in g/g
      &  xm_liq(klev),             & !< specific liquid water content
      &  xm_ice(klev),             & !< specific ice content in g/g
      &  cdnc(klev),               & !< cloud nuclei concentration
      &  cld_frc(klev),            & !< fractional cloud cover
      &  xm_o3(klev),              & !< o3 mass mixing ratio
      ! >> sylvia_20200305
      ! Aerosol inputs from external data.
      &  aotss(12),                & !< monthly climatology, sea salt AOT
      &  aotorg(12),               & !< monthly climatology, organic AOT
      &  aotbc(12),                & !< monthly climatology, black carbon AOT
      &  aotso4(12),               & !< monthly climatology, sulfate AOT
      &  aotdu(12)                   !< monthly climatology, dust AOT
      ! << sylvia_20200305

    REAL(wp), INTENT(out) ::        &
      &  flx_lw_net(klev+1),        & !< net downward LW flux profile,
      &  flx_sw_net(klev+1),        & !< net downward SW flux profile,
      &  flx_lw_net_clr(klev+1),    & !< clrsky downward LW flux profile,
      &  flx_sw_net_clr(klev+1),    & !< clrsky downward SW flux profile,
      &  flx_uplw_sfc,              & !< sfc LW upward flux,
      &  flx_upsw_sfc,              & !< sfc SW upward flux,
      &  flx_uplw_sfc_clr,          & !< clrsky sfc LW upward flux,
      &  flx_upsw_sfc_clr             !< clrsky sfc SW upward flux,

    REAL(wp), INTENT(out), OPTIONAL ::    &
      &  flx_upsw_toa                !< TOA SW upward flux,

    ! >> sylvia_20200302
    ! Set the number of shortwave and longwave spectral bands for RRTM.
    ! SW value can be verified within mo_srtm_config module as jpsw variable.
    ! LW value can be verified within mo_lrtm_par module as nbndlw variable.
    INTEGER, PARAMETER   :: jpsw   = 14
    INTEGER, PARAMETER   :: jpband = 16
    ! << sylvia_20200302
    
    ! >> sylvia_20200305
    ! Pulling parameters from shared/mo_impl_constants.f90
    ! identifiers for aerosol classes of Tegen climatology 
    INTEGER, PARAMETER :: iss   =  1
    INTEGER, PARAMETER :: iorg  =  2
    INTEGER, PARAMETER :: ibc   =  3
    INTEGER, PARAMETER :: iso4  =  4
    INTEGER, PARAMETER :: idu   =  5
    INTEGER, PARAMETER :: nclass_aero = 5
    ! << sylvia_20200305
    
    ! >> sylvia_20200303
    ! Pulling parameters from the mo_nwp_rrtm_interface module.
    REAL(wp), PARAMETER::  &
      & zaeops = 0.05_wp,   &
      & zaeopl = 0.2_wp,    &
      & zaeopu = 0.1_wp,    &
      & zaeopd = 1.9_wp,    &
      & ztrpt  = 30.0_wp,   &
      & ztrbga = 0.03_wp  / (101325.0_wp - 19330.0_wp), &
      & zvobga = 0.007_wp /  19330.0_wp , &
      & zstbga = 0.045_wp /  19330.0_wp
    ! << sylvia_20200303
    
    ! >> sylvia_20200303
    ! Pulling parameters from the aerdis subroutine in the mo_radiation_rg_par module
    REAL (wp), PARAMETER  ::  &
      zhss = 8434.0_wp/1000.0_wp ,  & !
      zhsl = 8434.0_wp/1000.0_wp ,  & !
      zhsu = 8434.0_wp/1000.0_wp ,  & !
      zhsd = 8434.0_wp/3000.0_wp      !
    REAL(wp) :: log_eta 
    ! << sylvia_20200303
    
    ! >> sylvia_20200305
    ! Changing imo1 and imo2 to constants for this simulation (previous and current months)
    ! zw is a constant weighting calculated from the seconds elapsed in the month.
    ! See calculate_time_interpolation_weights subroutine in the mo_bcs_time_interpolation module.
    INTEGER, PARAMETER :: imo1  =  7
    INTEGER, PARAMETER :: imo2  =  8
    REAL, PARAMETER    :: zw    =  0.758
    !for Tegen aerosol time interpolation
    ! << sylvia_20200305
    
    INTEGER  :: jk, jl, jp, jkb, jspec,   & !< loop indicies
      &  icldlyr(klev)                      !< index for clear or cloudy

    REAL(wp) ::                       &
      &  zsemiss(jpband),             & !< LW surface emissivity by band
      &  ppd_hl(klev),                & !< pressure thickness in Pa
      &  pm_sfc,                      & !< surface pressure in hPa
      &  amm,                         & !< molecular weight of moist air
      &  delta,                       & !< pressure thickness
      &  zscratch,                    & !< scratch array
      ! >> sylvia_20200305, zaeq* parameters added back in here.
      ! Aerosol array to hold what was prm_diag%aerosol.
      &  zaeq1(klev),                 & !< aerosol continental
      &  zaeq2(klev),                 & !< aerosol maritime
      &  zaeq3(klev),                 & !< aerosol urban
      &  zaeq4(klev),                 & !< aerosol volcano ashes
      &  zaeq5(klev),                 & !< aerosol stratospheric background
      &  aerosol(5)                     !< 5 classes of aerosol
      ! << sylvia_20200305

    REAL(wp) :: z_sum_aea, z_sum_aes !help variables for aerosol
    
    ! >> sylvia_20200305
    ! Pulling the intermediate variables for treatment of aerosol from nwp_ozon_aerosol.
    REAL(wp)::                        &
      & zsign(klev),                  &  ! looking at pvdaes in aerdis subroutine
      & zvdaes(klev),                 &
      & zvdael(klev),                 &
      & zvdaeu(klev),                 &
      & zvdaed(klev),                 &
      & zaetr_top, zaetr_bot, zaetr,  &
      & zaeqdo,    zaeqdn,            &
      & zaequo,    zaequn,            &
      & zaeqlo,    zaeqln,            &
      & zaeqso,    zaeqsn,    zw    ! (nproma,pt_patch%nblks_c)
    ! << sylvia_20200305
    
    !
    ! --- vertically reversed _vr variables
    !
    REAL(wp) ::                           &
      &  col_dry_vr(klev),          & !< number of molecules/cm2 of
      &  pm_fl_vr(klev),            & !< full level pressure [hPa]
      &  pm_hl_vr(klev+1),          & !< half level pressure [hPa]
      &  tk_fl_vr(klev),            & !< full level temperature [K]
      &  tk_hl_vr(klev+1),          & !< half level temperature [K]
      &  cdnc_vr(klev),             & !< cloud nuclei concentration
      &  cld_frc_vr(klev),          & !< secure cloud fraction
      &  ziwgkg_vr(klev),           & !< specific ice water content
      &  ziwc_vr(klev),             & !< ice water content per volume
      &  ziwp_vr(klev),             & !< ice water path in g/m2
      &  zlwgkg_vr(klev),           & !< specific liquid water content
      &  zlwp_vr(klev),             & !< liquid water path in g/m2
      &  zlwc_vr(klev),             & !< liquid water content per
      &  wkl_vr(jpinpx,klev),       & !< number of molecules/cm2 of
      &  wx_vr(jpxsec,klev),        & !< number of molecules/cm2 of
      &  cld_tau_lw_vr(klev,jpband),& !< LW optical thickness of clouds
      &  cld_tau_sw_vr(jpsw,klev),  & !< extincion
      &  cld_cg_sw_vr(jpsw,klev),   & !< asymmetry factor
      &  cld_piz_sw_vr(jpsw,klev),  & !< single scattering albedo
      &  aer_tau_lw_vr(klev,jpband),& !< LW optical thickness of aerosols
      &  aer_tau_sw_vr(klev,jpsw),  & !< aerosol optical thickness
      &  aer_cg_sw_vr(klev,jpsw),   & !< aerosol asymmetry factor
      &  aer_piz_sw_vr(klev,jpsw),  & !< aerosol single scattering albedo
      &  flx_uplw_vr(klev+1),       & !< upward flux, total sky
      &  flx_uplw_clr_vr(klev+1),   & !< upward flux, clear sky
      &  flx_dnlw_vr(klev+1),       & !< downward flux, total sky
      &  flx_dnlw_clr_vr(klev+1),   & !< downward flux, clear sky
      &  flx_upsw(klev+1),          & !< upward flux total sky
      &  flx_upsw_clr(klev+1),      & !< upward flux clear sky
      &  flx_dnsw(klev+1),          & !< downward flux total sky
      &  flx_dnsw_clr(klev+1)         !< downward flux clear sky
    ! >> sylvia_20200303, Removing the tune_dust parameter here.

    REAL(wp) ::                  &
         re_drop   (klev),       & !< effective radius of liquid
         re_cryst  (klev),       & !< effective radius of ice
         aux_out   (9),          &
         zmu0         ,          &
         zdayfrc

    INTEGER, PARAMETER    :: rng_seed_size = 4
    INTEGER :: rnseeds(rng_seed_size)
    
    ! >> sylvia_20200302
    ! Fix volume mixing ratios of trace gases according to namelist values.
    ! Below am* variables are taken from shared/mo_physical_constants.f90
    REAL(wp), PARAMETER   :: amd       = 28.970_wp        !> [g/mol] dry air
    REAL(wp), PARAMETER   :: amco2     = 44.011_wp        !>[g/mol] CO2
    REAL(wp), PARAMETER   :: amch4     = 16.043_wp        !! [g/mol] CH4
    REAL(wp), PARAMETER   :: amo3      = 47.9982_wp       !! [g/mol] O3
    REAL(wp), PARAMETER   :: amo2      = 31.9988_wp       !! [g/mol] O2
    REAL(wp), PARAMETER   :: amn2o     = 44.013_wp        !! [g/mol] N2O
    REAL(wp), PARAMETER   :: amc11     = 137.3686_wp      !! [g/mol] CFC11
    REAL(wp), PARAMETER   :: amc12     = 120.9140_wp      !! [g/mol] CFC12
    REAL(wp), PARAMETER   :: amw       = 18.0154_wp       !! [g/mol] H2O
    REAL, PARAMETER       :: vmr_co2   = 390.e-6_wp       !! [ppmv] CO2
    REAL, PARAMETER       :: vmr_ch4   = 1800.e-9_wp      !! [ppmv] CH4
    REAL, PARAMETER       :: vmr_n2o   = 322.0e-9_wp      !! [ppmv] N20
    REAL, PARAMETER       :: vmr_o2    = 0.20946_wp       !! [ppmv] O2
    REAL, PARAMETER       :: vmr_cfc11 = 240.e-12_wp      !! [ppmv] CFC11
    REAL, PARAMETER       :: vmr_cfc12 = 532.3e-12_wp     !! [ppmv] CFC12
    ! << sylvia_20200302
    
    ! >> sylvia_20200302
    ! Add in the variables necessary for the imported newcld_optics subroutine.
    INTEGER  :: iband, jk, jl, ml1, ml2, mi1, mi2
    REAL(wp) :: ztol, ztoi, zol, zoi, zgl, zgi, wl1, wl2, wi1, wi2, zfact
    REAL(wp) :: zmsald, zmsaid
    REAL(wp) :: &
    &  zkap,        & !< breath parameter for scaling effective radius
    &  zlwpt,       & !< liquid water path
    &  zinhoml,     & !< cloud inhomogeneity factor (liquid)
    &  re_droplets,        & !< effective radius of liquid distribution
    &  re_crystals,        & !< effective radius of ice distribution
    &  zscratch,           & !< effective radius of ice distribution
    &  ztau(klev,n_mdl_bnds), & !< SW optical depth
    &  zomg(klev,n_mdl_bnds), & !< cloud single scattering albedo
    &  zasy(klev,n_mdl_bnds)    !< cloud asymmetry factor
    ! << sylvia_20200302
    
    ! Initialize output variables
    flx_lw_net(:,:)     = 0._wp
    flx_lw_net_clr(:,:) = 0._wp
    flx_sw_net(:,:)     = 0._wp
    flx_sw_net_clr(:,:) = 0._wp
    flx_uplw_sfc(:)     = 0._wp
    flx_uplw_sfc_clr(:) = 0._wp
    flx_upsw_sfc(:)     = 0._wp
    flx_upsw_sfc_clr(:) = 0._wp
    IF (PRESENT(flx_dnsw_diff_sfc)) flx_dnsw_diff_sfc(:) = 0._wp
    IF (PRESENT(flx_upsw_toa))      flx_upsw_toa(:)      = 0._wp
    IF (PRESENT(flx_dnpar_sfc))     flx_dnpar_sfc(:)     = 0._wp
    IF (PRESENT(vis_frc_sfc))       vis_frc_sfc(:)       = 0._wp
    IF (PRESENT(nir_dff_frc_sfc))   nir_dff_frc_sfc(:)   = 0._wp
    IF (PRESENT(vis_dff_frc_sfc))   vis_dff_frc_sfc(:)   = 0._wp
    IF (PRESENT(par_dff_frc_sfc))   par_dff_frc_sfc(:)   = 0._wp

    ! >> sylvia_20200302
    ! Calculate the mass mixing ratios from volume mixing ratios.
    xm_co2(klev)      = vmr_o2*amco2/amd         !< [g/g] co2 mass mixing ratio
    xm_ch4(klev)      = vmr_ch4*amch4/amd        !< [g/g] ch4 mass mixing ratio
    xm_n2o(klev)      = vmr_n2o*amn2o/amd        !< [g/g] n2o mass mixing ratio
    xm_cfc11(klev)    = vmr_cfc11*amc11/amd      !< [g/g] cfc11 mass mixing ratio
    xm_cfc12(klev)    = vmr_cfc12*amc12/amd      !< [g/g] cfc 12 mass mixing ratio
    xm_o2(klev)       = vmr_o2*amo2/amd          !< [g/g] o2  mass mixing ratio
    ! << sylvia_20200302
    
    
    ! 1.0 Constituent properties
    !--------------------------------
    ! Control for infintesimal cloud fractions
    DO jk = 1, klev
      jkb = klev+1-jk
      cld_frc_vr(jk)  = cld_frc(jkb)

      IF (cld_frc_vr(jk) > 2.0_wp*EPSILON(1.0_wp)) THEN
        ! only clouds > 2 epsilon are made visible to radiation
        icldlyr  (jk) = 1
        ziwgkg_vr(jk) = xm_ice(jkb)*1000.0_wp/cld_frc_vr(jk)
        zlwgkg_vr(jk) = xm_liq(jkb)*1000.0_wp/cld_frc_vr(jk)
      ELSE
        ! clouds <= 2 epsilon are ade invisble to radiation
        icldlyr  (jk) = 0
        ziwgkg_vr(jk) = 0.0_wp
        zlwgkg_vr(jk) = 0.0_wp
      ENDIF
    END DO
    !
    ! --- main constituent reordering
    !
    pm_hl_vr(klev+1) = 0.01_wp*pp_hl(1)
    tk_hl_vr(klev+1) = tk_hl(1)
    pm_sfc          = 0.01_wp*pp_sfc

    DO jk = 1, klev
      jkb = klev+1-jk
      ! initialization
      wkl_vr(:,:,jk) = 0.0_wp
      wx_vr (:,:,jk) = 0.0_wp
      delta = pp_hl(jkb+1)-pp_hl(jkb)
      !
      ! --- thermodynamic arrays
      !
      pm_hl_vr(jk) = 0.01_wp*pp_hl(jkb+1)
      pm_fl_vr(jk) = 0.01_wp*pp_fl(jkb)
      tk_hl_vr(jk) = tk_hl(jkb+1)
      tk_fl_vr(jk) = tk_fl(jkb)
      !
      ! --- cloud properties
      !
      zscratch       = pp_fl(jkb)/tk_fl(jkb)
      ziwc_vr(jk) = ziwgkg_vr(jk)*zscratch/rd
      ziwp_vr(jk) = ziwgkg_vr(jk)*delta/grav
      zlwc_vr(jk) = zlwgkg_vr(jk)*zscratch/rd
      zlwp_vr(jk) = zlwgkg_vr(jk)*delta/grav
      cdnc_vr(jk) = cdnc(jkb)*1.e-6_wp
      !
      ! --- radiatively active gases
      !
      wkl_vr(1,jk)   = xm_vap(jkb)*amd/amw
      wkl_vr(2,jk)   = xm_co2(jkb)*amd/amco2
      wkl_vr(3,jk)   = xm_o3(jkb) *amd/amo3
      wkl_vr(4,jk)   = xm_n2o(jkb)*amd/amn2o
      wkl_vr(6,jk)   = xm_ch4(jkb)*amd/amch4
      wkl_vr(7,jk)   = xm_o2 (jkb)*amd/amo2
      amm               = (1.0_wp-wkl_vr(1,jk))*amd + wkl_vr(1,jk)*amw
      col_dry_vr(jk) = (0.01_wp*delta)*10.0_wp*avo/grav/amm / (1.0_wp+wkl_vr(1,jk))
      !
      ! --- alternate treatment for cfcs
      !
      wx_vr(2,jk) = col_dry_vr(jk)*xm_cfc11(jkb)*1.e-20_wp
      wx_vr(3,jk) = col_dry_vr(jk)*xm_cfc12(jkb)*1.e-20_wp
  END DO
  
  DO jp = 1, 7
    wkl_vr(jp,:)=col_dry_vr(:)*wkl_vr(jp,:)
  END DO
  !
  ! 2.0 Surface Properties
  ! --------------------------------
  zsemiss(:) = SPREAD(emis_rad,2,jpband)
  !
  ! 3.0 Particulate Optical Properties
  ! --------------------------------
  ppd_hl(:) = pp_hl(2:klev+1)-pp_hl(1:klev)

  !
  ! >> sylvia_20200305 START: nwp_ozon_aerosol
  !
  ! >> sylvia_20200305
  ! Changing prm_diag%aerosol nomenclautre to just aerosol array with constant
  ! indices and weightings.
  aerosol(iss)  = aotss(imo1) + ( aotss(imo2) - aotss(imo1)  ) * zw
  aerosol(iorg) = aotorg(imo1) + ( aotorg(imo2)  - aotorg(imo1)  ) * zw
  aerosol(ibc)  = aotbc(imo1) + ( aotbc(imo2)  -  aotbc(imo1)   ) * zw
  aerosol(iso4) = aotso4(imo1) +  ( aotso4(imo2)  - aotso4(imo1)  ) * zw
  aerosol(idu)  = aotdu(imo1) + ( aotdu(imo2) - aotdu(imo1) ) * zw
  ! << sylvia_20200305
  
  ! Below changed from pt_diag%pres_ifc(jk,jb) to pp_hl
  DO jk = 2, klev
     zsign(jk) = pp_hl(jk) / 101325._wp 
  ENDDO
  
  ! The routine aerdis is called to receive some parameters for the vertical
  ! distribution of background aerosol.
  ! >> sylvia_20200304
  ! Adding the aerdis subroutine code in here.
  zvdaes(1) = 0.0_wp
  zvdael(1) = 0.0_wp
  zvdaeu(1) = 0.0_wp
  zvdaed(1) = 0.0_wp
    
  DO jk = 2,klev
     log_eta    = LOG(petah(jk))
     zvdaes(jk) = EXP(zhss*log_eta) ! petah(jc,jk)**zhss
     zvdael(jk) = zsign(jk)     ! petah(jc,jk)**zhsl; zhsl is the same as zhss
     zvdaeu(jk) = zsign(jk)     ! petah(jc,jk)**zhsu; zhsu is the same as zhss
     zvdaed(jk) = EXP(zhsd*log_eta) ! petah(jc,jk)**zhsd
  ENDDO
  ! << sylvia_20200304
  
  ! loop over layers
  DO jk = 1,klev
     zaeqsn       = aerosol(iss)                     * zvdaes(jk+1)
     zaeqln       = (aerosol(iorg) + aerosol(iso4))  * zvdael(jk+1)
     zaequn       = aerosol(ibc)                     * zvdaeu(jk+1)
     zaeqdn       = aerosol(idu)                     * zvdaed(jk+1)
     zaetr_bot    = zaetr_top &
            & * ( MIN (1.0_wp, tk_hl(jk)/tk_hl(jk+1)) )**ztrpt

     zaetr        = SQRT(zaetr_bot*zaetr_top)
     zaeq1(jk)    = (1.0_wp - zaetr)*( ztrbga*pt_diag%dpres_mc(jk) &
            &            + zaeqln - zaeqlo )
     zaeq2(jk)    = (1.0_wp - zaetr) * (zaeqsn - zaeqso)
     zaeq3(jk)    = (1.0_wp - zaetr) * (zaeqdn - zaeqdo)
     zaeq4(jk)    = (1.0_wp - zaetr) * (zaequn - zaequo)
     zaeq5(jk)    = zaetr  *  zstbga  *  pt_diag%dpres_mc(jk)

     zaetr_top    = zaetr_bot
     zaeqso       = zaeqsn
     zaeqlo       = zaeqln
     zaequo       = zaequn
     zaeqdo       = zaeqdn
  ENDDO
  ! << sylvia_20200305 STOP: nwp_ozon_aerosol
  !
  ! >> sylvia_20200303
  ! Adding in variables for zaea_rrtm (absorption), zaes_rrtm (scattering), and 
  ! zaeg_rrtm (asymmetry factor). zaea_rrtm is the "ratio of optical thickness
  ! for the absorption in spectral interval jpspec and total optical thickness
  ! at 0.55 um for an aerosoltyp specified by the second array index. zaes_rrtm
  ! is the "analog for the optical thickness of scattering." These values are 
  ! pulled from mo_aerosol_util SUBROUTINE init_aerosol_props_tegen_rrtm.
  
  ! the following aerosol types (second array index) are considered:
  ! 1. continental, 2. maritime, 3. desert, 4. urban, 5. stratospheric background (SB)

   !absorption
   zaea_rrtm=RESHAPE( (/ &
     &0.0304_wp,0.0367_wp,0.0462_wp,0.0566_wp,0.0496_wp,0.0336_wp,0.0355_wp,0.0456_wp,&
     &0.0272_wp,0.0264_wp,0.0290_wp,0.0156_wp,0.0165_wp,0.0157_wp,0.0138_wp,0.0401_wp,&
     &0.0401_wp,0.0760_wp,0.0214_wp,0.0227_wp,0.0295_wp,0.0394_wp,0.0431_wp,0.0519_wp,&
     &0.0611_wp,0.0774_wp,0.1012_wp,0.1412_wp,0.2632_wp,0.0324_wp,                    &
     &0.1096_wp,0.1614_wp,0.2294_wp,0.2506_wp,0.2242_wp,0.1190_wp,0.0680_wp,0.0664_wp,&
     &0.0656_wp,0.0749_wp,0.1250_wp,0.0425_wp,0.0498_wp,0.0425_wp,0.0259_wp,0.1619_wp,&
     &0.1619_wp,0.2152_wp,0.0139_wp,0.0119_wp,0.0046_wp,0.0036_wp,0.0020_wp,0.0016_wp,&
     &0.0012_wp,0.0013_wp,0.0016_wp,0.0035_wp,0.0147_wp,0.0882_wp,                    &
     &0.0974_wp,0.1529_wp,0.1643_wp,0.1373_wp,0.1753_wp,0.1923_wp,0.2804_wp,0.2426_wp,&
     &0.1263_wp,0.1321_wp,0.0979_wp,0.0664_wp,0.0360_wp,0.0311_wp,0.0325_wp,0.0833_wp,&
     &0.0833_wp,0.1170_wp,0.0739_wp,0.0631_wp,0.0604_wp,0.0628_wp,0.0645_wp,0.0677_wp,&
     &0.0843_wp,0.1328_wp,0.2224_wp,0.3022_wp,0.3579_wp,0.1820_wp,                    &
     &0.0267_wp,0.0329_wp,0.0420_wp,0.0515_wp,0.0461_wp,0.0332_wp,0.0354_wp,0.0447_wp,&
     &0.0303_wp,0.0306_wp,0.0342_wp,0.0248_wp,0.0274_wp,0.0276_wp,0.0271_wp,0.0526_wp,&
     &0.0526_wp,0.0903_wp,0.0450_wp,0.0492_wp,0.0596_wp,0.0754_wp,0.0842_wp,0.1082_wp,&
     &0.1429_wp,0.1926_wp,0.2595_wp,0.3379_wp,0.4761_wp,0.0340_wp,                    &
     &0.0060_wp,0.0117_wp,0.0269_wp,0.0222_wp,0.0195_wp,0.0398_wp,0.0733_wp,0.1091_wp,&     ! SB
     &0.1124_wp,0.0415_wp,0.0424_wp,0.0495_wp,0.0451_wp,0.0484_wp,0.0540_wp,0.0735_wp,&     ! SB
     &0.0735_wp,0.0188_wp,0.0021_wp,0.0014_wp,0.0007_wp,0.0002_wp,0.0000_wp,0.0000_wp,&     ! SB
     &0.0000_wp,0.0000_wp,0.0000_wp,0.0000_wp,0.0000_wp,0.0628_wp/),(/jpsw+jpband,5/))      ! SB
     
  !scattering
   zaes_rrtm=RESHAPE( (/ &
     &0.0060_wp,0.0107_wp,0.0134_wp,0.0150_wp,0.0152_wp,0.0200_wp,0.0232_wp,0.0211_wp,&
     &0.0112_wp,0.0186_wp,0.0128_wp,0.0260_wp,0.0339_wp,0.0368_wp,0.0409_wp,0.0527_wp,&
     &0.0527_wp,0.0621_wp,0.0715_wp,0.0929_wp,0.1276_wp,0.1895_wp,0.2350_wp,0.3930_wp,&
     &0.6641_wp,0.9834_wp,1.3737_wp,1.7160_wp,1.9115_wp,0.0198_wp,                    &
     &0.0188_wp,0.0421_wp,0.0576_wp,0.0547_wp,0.0430_wp,0.0367_wp,0.0806_wp,0.1209_wp,&
     &0.1681_wp,0.2257_wp,0.2440_wp,0.3622_wp,0.4540_wp,0.5026_wp,0.5765_wp,0.5986_wp,&
     &0.5986_wp,0.5225_wp,0.7420_wp,0.8311_wp,0.8970_wp,0.9444_wp,0.9637_wp,0.9763_wp,&
     &0.9855_wp,1.0034_wp,1.0337_wp,1.0640_wp,1.0795_wp,0.1312_wp,                    &
     &0.0458_wp,0.0823_wp,0.0667_wp,0.0642_wp,0.1080_wp,0.1471_wp,0.2422_wp,0.1216_wp,&
     &0.0717_wp,0.1616_wp,0.2027_wp,0.3042_wp,0.4045_wp,0.4369_wp,0.4685_wp,0.5043_wp,&
     &0.5043_wp,0.5782_wp,0.6898_wp,0.7477_wp,0.7926_wp,0.8320_wp,0.8503_wp,0.8736_wp,&
     &0.8874_wp,0.8737_wp,0.8278_wp,0.7857_wp,0.7571_wp,0.1714_wp,                    &
     &0.0048_wp,0.0085_wp,0.0107_wp,0.0119_wp,0.0121_wp,0.0160_wp,0.0185_wp,0.0170_wp,&
     &0.0090_wp,0.0150_wp,0.0103_wp,0.0210_wp,0.0274_wp,0.0298_wp,0.0332_wp,0.0430_wp,&
     &0.0430_wp,0.0485_wp,0.0593_wp,0.0776_wp,0.1073_wp,0.1610_wp,0.2008_wp,0.3398_wp,&
     &0.5809_wp,0.8701_wp,1.2309_wp,1.5535_wp,1.7368_wp,0.0159_wp,                    &
     &0.0000_wp,0.0000_wp,0.0000_wp,0.0000_wp,0.0001_wp,0.0003_wp,0.0006_wp,0.0008_wp,&     ! SB
     &0.0005_wp,0.0003_wp,0.0008_wp,0.0013_wp,0.0024_wp,0.0030_wp,0.0040_wp,0.0059_wp,&     ! SB
     &0.0059_wp,0.0123_wp,0.0236_wp,0.0384_wp,0.0651_wp,0.1246_wp,0.1801_wp,0.3807_wp,&     ! SB
     &0.7105_wp,1.0514_wp,1.3754_wp,1.5334_wp,1.5495_wp,0.0009_wp/),(/jpsw+jpband,5/))      ! SB

   !asymmetry factor
   zaeg_rrtm=RESHAPE( (/ &
     &0.4388_wp,0.5396_wp,0.6191_wp,0.6535_wp,0.6876_wp,0.6718_wp,0.6493_wp,0.6782_wp,&
     &0.7958_wp,0.7537_wp,0.7757_wp,0.7821_wp,0.7583_wp,0.7487_wp,0.7351_wp,0.6917_wp,&
     &0.6917_wp,0.6989_wp,0.6982_wp,0.6726_wp,0.6426_wp,0.6294_wp,0.6337_wp,0.6582_wp,&
     &0.6850_wp,0.7061_wp,0.7212_wp,0.7306_wp,0.7417_wp,0.6978_wp,                    &
     &0.4062_wp,0.4507_wp,0.4878_wp,0.5302_wp,0.5850_wp,0.6962_wp,0.7242_wp,0.7293_wp,&
     &0.7414_wp,0.7484_wp,0.7607_wp,0.7785_wp,0.7805_wp,0.7785_wp,0.7724_wp,0.7690_wp,&
     &0.7690_wp,0.8348_wp,0.8316_wp,0.8170_wp,0.8074_wp,0.7990_wp,0.7954_wp,0.7897_wp,&
     &0.7884_wp,0.7927_wp,0.8001_wp,0.8057_wp,0.8076_wp,0.7462_wp,                    &
     &0.4219_wp,0.3928_wp,0.5306_wp,0.6229_wp,0.5544_wp,0.5454_wp,0.4353_wp,0.5736_wp,&
     &0.7502_wp,0.6957_wp,0.7038_wp,0.6881_wp,0.6740_wp,0.6739_wp,0.6784_wp,0.6969_wp,&
     &0.6969_wp,0.7068_wp,0.6965_wp,0.6918_wp,0.6904_wp,0.6911_wp,0.6915_wp,0.6952_wp,&
     &0.7080_wp,0.7326_wp,0.7689_wp,0.8000_wp,0.8206_wp,0.5788_wp,                    &
     &0.4387_wp,0.5394_wp,0.6187_wp,0.6531_wp,0.6871_wp,0.6712_wp,0.6482_wp,0.6756_wp,&
     &0.7930_wp,0.7498_wp,0.7685_wp,0.7766_wp,0.7520_wp,0.7419_wp,0.7277_wp,0.6828_wp,&
     &0.6828_wp,0.6875_wp,0.6872_wp,0.6622_wp,0.6333_wp,0.6209_wp,0.6250_wp,0.6479_wp,&
     &0.6725_wp,0.6912_wp,0.7043_wp,0.7129_wp,0.7254_wp,0.6956_wp,                    &
     &0.0021_wp,0.0039_wp,0.0061_wp,0.0078_wp,0.0109_wp,0.0161_wp,0.0201_wp,0.0206_wp,&     ! SB
     &0.0217_wp,0.0320_wp,0.0428_wp,0.0583_wp,0.0773_wp,0.0856_wp,0.0985_wp,0.1310_wp,&     ! SB
     &0.1310_wp,0.1906_wp,0.2625_wp,0.3154_wp,0.3869_wp,0.4787_wp,0.5279_wp,0.6272_wp,&     ! SB
     &0.6941_wp,0.7286_wp,0.7358_wp,0.7177_wp,0.6955_wp,0.0616_wp/),(/jpsw+jpband,5/))      ! SB   
  ! << sylvia_20200303
    
  ! >> sylvia_20200303 
  ! Below I have taken out tune_dust for simplicity, assuming it equals one.
  ! It would multiply the third term.
  DO jspec=1,jpband
    DO jk=1,klev
      jkb = klev+1-jk
      ! LW opt thickness of aerosols
      aer_tau_lw_vr(jk,jspec) =  zaeq1(jkb) * zaea_rrtm(jspec,1) &
      &                         + zaeq2(jl,jkb) * zaea_rrtm(jspec,2) &
      &                         + zaeq3(jl,jkb) * zaea_rrtm(jspec,3) &
      &                         + zaeq4(jl,jkb) * zaea_rrtm(jspec,4) &
      &                         + zaeq5(jl,jkb) * zaea_rrtm(jspec,5)
      
    ENDDO
  ENDDO
  
  DO jspec=1+jpband,jpband+jpsw
    DO jk=1,klev
      jkb = klev+1-jk
      
      z_sum_aea = zaeq1(jkb) * zaea_rrtm(jspec,1) &
      &       + zaeq2(jkb) * zaea_rrtm(jspec,2) &
      &       + zaeq3(jkb) * zaea_rrtm(jspec,3) &
      &       + zaeq4(jkb) * zaea_rrtm(jspec,4) &
      &       + zaeq5(jkb) * zaea_rrtm(jspec,5)

      z_sum_aes = zaeq1(jkb) * zaes_rrtm(jspec,1) &
      &       + zaeq2(jkb) * zaes_rrtm(jspec,2) &
      &       + zaeq3(jkb) * zaes_rrtm(jspec,3) &
      &       + zaeq4(jkb) * zaes_rrtm(jspec,4) &
      &       + zaeq5(jkb) * zaes_rrtm(jspec,5)

      ! sw aerosol optical thickness
      aer_tau_sw_vr(jl,jk,jspec-jpband) = z_sum_aea + z_sum_aes

      ! sw aerosol single scattering albedo
      aer_piz_sw_vr(jl,jk,jspec-jpband) = z_sum_aes / ( z_sum_aea + z_sum_aes )

      ! sw aerosol asymmetry factor
      aer_cg_sw_vr(jl,jk,jspec-jpband) =                                  &
      & (   zaeq1(jl,jkb) * zaes_rrtm(jspec,1) * zaeg_rrtm(jspec,1)   &
      &   + zaeq2(jl,jkb) * zaes_rrtm(jspec,2) * zaeg_rrtm(jspec,2)   &
      &   + zaeq3(jl,jkb) * zaes_rrtm(jspec,3) * zaeg_rrtm(jspec,3)   &
      &   + zaeq4(jl,jkb) * zaes_rrtm(jspec,4) * zaeg_rrtm(jspec,4)   &
      &   + zaeq5(jl,jkb) * zaes_rrtm(jspec,5) * zaeg_rrtm(jspec,5) ) / z_sum_aes
    ENDDO
  ENDDO

  ! >> sylvia_20200302
  ! Replacing the call to newcld_optics.f90 with the code therein.
    ! First check for consistency of number of LW and SW bands
  IF (nbnds_lw /= jpband .OR. nbnds_sw /= jpsw) THEN
    CALL finish('mo_newcld_optics: Inconsistent number of spectral bands')
  END IF

  IF (relmin < 1.5_wp .OR. relmin > 2.5_wp ) THEN
    CALL finish('Apparently unsuccessful loading of optical tables')
  END IF
  !
  ! 1.0 Basic cloud properties
  ! --------------------------------
  zinhoml = 0.77_wp
  zkap    = zkap_cont*(zland) + zkap_mrtm*(1.0_wp-zland)
  !
  ! 2.0 Cloud Optical Properties by interpolating tables in effective radius
  ! --------------------------------

  zfact = 1.0e6_wp*(3.0e-9_wp/(4.0_wp*pi*rhoh2o))**(1.0_wp/3.0_wp)
  DO jk=1,klev
      IF (icldlyr(jk)==1 .AND. (zlwp_vr(jk)+ziwp_vr(jk))>ccwmin) THEN

        ! see ECHAM5 documentation (Roeckner et al, MPI report 349)
        re_crystals = MAX(reimin ,MIN(reimax  ,83.8_wp*ziwc_vr(jk)**0.216_wp))
        re_droplets = MAX(relmin,MIN(relmax,zfact*zkap(jl)*(zlwc_vr(jk) &
          & /cdnc_vr(jk))**(1.0_wp/3.0_wp)))

        ! alternative formulation Ou and Liou (1995) as function of T as in IFS (cy38r2)
        ! re_crystals = MAX(20.0_wp, MIN(70.0_wp, 0.5_wp * &  ! limits to range of data
        !   & ( 326.3_wp                                   &
        !   & + 12.42_wp  * (tk_fl_vr(jk)-tmelt)              &
        !   & + 0.197_wp  * (tk_fl_vr(jk)-tmelt)**2           &
        !   & + 0.0012_wp * (tk_fl_vr(jk)-tmelt)**3 )))

        ! optional tuning of effective crystal radius
        ! re_crystals = re_crystals * 1.25_wp

        ml1 = MAX(1,MIN(n_sizes-1,FLOOR(1.0_wp+(re_droplets-relmin)/del_rel)))
        ml2 = ml1 + 1
        wl1 = 1.0_wp - (re_droplets - (relmin + del_rel* REAL(ml1-1,wp)) )/del_rel
        wl2 = 1.0_wp - wl1

        mi1 = MAX(1,MIN(n_sizes-1,FLOOR(1.0_wp+(re_crystals-reimin)/del_rei)))
        mi2 = mi1 + 1
        wi1 = 1.0_wp - (re_crystals - (reimin + del_rei * REAL(mi1-1,wp)) )/del_rei
        wi2 = 1.0_wp - wi1

        ! Note: The following implementation is in principle intended; the first SW
        !       band overlaps with the last LW band, and the last SW band is unused
        !  ***  An open question is, however, is the usage of the coefficients
        !       in the overlapping band: currently, ztau and zomg are set for LW,
        !       whereas zasy are set for SW. ***
!cdir expand=nbnds_sw
        DO iband = nbnds_lw,n_mdl_bnds-1
          ztol = zlwp_vr(jk)*(wl1*z_ext_l(ml1,iband) + wl2*z_ext_l(ml2,iband))
          ztoi = ziwp_vr(jk)*(wi1*z_ext_i(mi1,iband) + wi2*z_ext_i(mi2,iband))
          zol  = 1.0_wp - (wl1*z_coa_l(ml1,iband) + wl2*z_coa_l(ml2,iband))
          zoi  = 1.0_wp - (wi1*z_coa_i(mi1,iband) + wi2*z_coa_i(mi2,iband))
          zgl  = wl1*z_asy_l(ml1,iband) + wl2*z_asy_l(ml2,iband)
          zgi  = wi1*z_asy_i(mi1,iband) + wi2*z_asy_i(mi2,iband)

          zscratch = (ztol*zol+ztoi*zoi)
          ztau(jl,jk,iband) = ztol*zinhoml(jl) + ztoi*zinhomi
          zomg(jl,jk,iband) = zscratch/(ztol+ztoi)
          zasy(jl,jk,iband) = (ztol*zol*zgl+ztoi*zoi*zgi)/zscratch
        END DO
        !
        ! set old Cloud Optics instead of Kinne Optics for LW
        !
!cdir expand=nbnds_lw
        DO iband = 1,nbnds_lw
          zmsald=0.025520637_wp+0.2854650784_wp*EXP(-0.088968393014_wp  &
               *re_droplets)
          zmsaid=(rebcuh(iband)+rebcug(iband)/re_crystals)
          ztau(jl,jk,iband)  = zmsald*zlwp_vr(jk)*zinhoml(jl)          &
               &             + zmsaid*ziwp_vr(jk)*zinhomi
          zomg(jl,jk,iband) = 0.0_wp
        END DO

        ELSE
!cdir begin expand=n_mdl_bnds-1
        ztau(jl,jk,1:n_mdl_bnds-1)  = 0.0_wp
        zomg(jl,jk,1:n_mdl_bnds-1)  = 1.0_wp
        zasy(jl,jk,1:n_mdl_bnds-1)  = 0.0_wp
!cdir end
      END IF
!cdir expand=nbnds_lw
      cld_tau_lw_vr(jl,jk,1:nbnds_lw) = ztau(jl,jk,1:nbnds_lw) * (1.0_wp - zomg(jl,jk,1:nbnds_lw))
!cdir begin expand=nbnds_sw
      cld_tau_sw_vr(jl,1:nbnds_sw,jk) = ztau(jl,jk,nbnds_lw:nbnds_lw+nbnds_sw-1)
      cld_piz_sw_vr(jl,1:nbnds_sw,jk) = zomg(jl,jk,nbnds_lw:nbnds_lw+nbnds_sw-1)
      cld_cg_sw_vr(jl,1:nbnds_sw,jk)  = zasy(jl,jk,nbnds_lw:nbnds_lw+nbnds_sw-1)
!cdir end
  END DO
  ! << sylvia_20200302

  !
  ! 4.0 Radiative Transfer Routines
  ! --------------------------------
  CALL lrtm(                                                                &
  !    input
  &    pm_fl_vr        ,pm_sfc          ,tk_fl_vr        ,tk_hl_vr       ,&
  &    tk_sfc          ,wkl_vr          ,wx_vr           ,col_dry_vr     ,&
  &    zsemiss         ,cld_frc_vr      ,cld_tau_lw_vr   ,aer_tau_lw_vr  ,&
  !    output
  &    flx_uplw_vr     ,flx_dnlw_vr     ,flx_uplw_clr_vr,flx_dnlw_clr_vr )


  CALL srtm_srtm_224gp(                                                     &
  !    input
  &    klev            ,jpsw                                             ,&
  &    alb_vis_dir     ,alb_nir_dir     ,alb_vis_dif     ,alb_nir_dif    ,&
  &    pm_fl_vr        ,tk_fl_vr        ,pmu0                            ,&
  &    col_dry_vr      ,wkl_vr                                           ,&
  &    cld_frc_vr      ,cld_tau_sw_vr   ,cld_cg_sw_vr    ,cld_piz_sw_vr  ,&
  &    aer_tau_sw_vr   ,aer_cg_sw_vr    ,aer_piz_sw_vr                   ,&
  &    ssi_radt                                                          ,&
  !    output
  &    flx_dnsw        ,flx_upsw        ,flx_dnsw_clr    ,flx_upsw_clr    )


  ! 5.0 Post Processing
  ! --------------------------------

  DO jk = 1, klev+1
    jkb = klev+2-jk
    flx_lw_net(jk)     = flx_dnlw_vr(jkb)-flx_uplw_vr(jkb)
    flx_lw_net_clr(jk) = flx_dnlw_clr_vr(jkb)-flx_uplw_clr_vr(jkb)
    flx_sw_net(jk)     = flx_dnsw(jk) - flx_upsw(jk)
    flx_sw_net_clr(jk) = flx_dnsw_clr(jk)-flx_upsw_clr(jk)
  END DO
    flx_uplw_sfc     = flx_uplw_vr(1)
    flx_uplw_sfc_clr = flx_uplw_clr_vr(1)
    flx_upsw_sfc     = flx_upsw(klev+1)
    flx_upsw_sfc_clr = flx_upsw_clr(klev+1)
    IF (PRESENT(flx_upsw_toa)) flx_upsw_toa = flx_upsw(1)

  END SUBROUTINE rrtm_interface
