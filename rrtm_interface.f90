  SUBROUTINE rrtm_interface(                                              &
    ! input                                                     ,&
    & jg              ,jb              ,klev                             ,&
    & ktype           ,zland           ,zglac                            ,&
    & pmu0                                                               ,&
    & alb_vis_dir     ,alb_nir_dir     ,alb_vis_dif     ,alb_nir_dif     ,&
    & emis_rad                                                           ,&
    & pp_fl           ,pp_hl           ,pp_sfc          ,tk_fl           ,&
    & tk_hl           ,tk_sfc          ,xm_vap                           ,&
    & xm_liq          ,xm_ice                                            ,&
    & cdnc                                                               ,&
    & cld_frc                                                            ,&
    & xm_o3           ,xm_co2          ,xm_ch4                           ,&
    & xm_n2o          ,xm_cfc11        ,xm_cfc12        ,xm_o2           ,&
    & zaeq1, zaeq2, zaeq3, zaeq4, zaeq5                                  ,&
    ! output
    & flx_lw_net      ,flx_sw_net      ,flx_lw_net_clr  ,flx_sw_net_clr  ,&
    & flx_uplw_sfc    ,flx_upsw_sfc    ,flx_uplw_sfc_clr,flx_upsw_sfc_clr,&
    ! optional output
    & flx_dnsw_diff_sfc, flx_upsw_toa  ,flx_dnpar_sfc                    ,&
    & vis_frc_sfc     ,nir_dff_frc_sfc ,vis_dff_frc_sfc ,par_dff_frc_sfc  )
    
    INTEGER,INTENT(in)  ::                &
      &  jg,                              & !< domain index
      &  jb,                              & !< block index
      &  klev                               !< number of levels

    INTEGER,INTENT(in)  ::                &
      &  ktype                       !< type of convection

    REAL(wp),INTENT(in) ::                &
      &  zland,                    & !< land-sea mask. (1. = land, 0. = sea/lakes)
      &  zglac,                    & !< fraction of land covered by glaciers
      &  pmu0,                     & !< mu0 for solar zenith angle
      &  alb_vis_dir,              & !< surface albedo for vis range and dir light
      &  alb_nir_dir,              & !< surface albedo for NIR range and dir light
      &  alb_vis_dif,              & !< surface albedo for vis range and dif light
      &  alb_nir_dif,              & !< surface albedo for NIR range and dif light
      &  emis_rad,                 & !< longwave surface emissivity
      &  pp_fl(klev),               & !< full level pressure in Pa
      &  pp_hl(klev+1),             & !< half level pressure in Pa
      &  pp_sfc(kbdim),                   & !< surface pressure in Pa
      &  tk_fl(klev),               & !< full level temperature in K
      &  tk_hl(klev+1),             & !< half level temperature in K
      &  tk_sfc(kbdim),                   & !< surface temperature in K
      &  xm_vap(klev),              & !< specific humidity in g/g
      &  xm_liq(klev),              & !< specific liquid water content
      &  xm_ice(klev),              & !< specific ice content in g/g
      &  cdnc(klev),                & !< cloud nuclei concentration
      &  cld_frc(klev),             & !< fractional cloud cover
      &  xm_o3(klev),               & !< o3 mass mixing ratio
      &  xm_co2(klev),              & !< co2 mass mixing ratio
      &  xm_ch4(klev),              & !< ch4 mass mixing ratio
      &  xm_n2o(klev),              & !< n2o mass mixing ratio
      &  xm_cfc11(klev),            & !< cfc 11 volume mixing ratio
      &  xm_cfc12(klev),            & !< cfc 12 volume mixing ratio
      &  xm_o2(klev),               & !< o2  mass mixing ratio
      &  zaeq1(klev),               & !< aerosol continental
      &  zaeq2(klev),               & !< aerosol maritime
      &  zaeq3(klev),               & !< aerosol urban
      &  zaeq4(klev),               & !< aerosol volcano ashes
      &  zaeq5(klev)                  !< aerosol stratospheric background


    REAL(wp), INTENT(out) ::              &
      &  flx_lw_net(klev+1),        & !< net downward LW flux profile,
      &  flx_sw_net(klev+1),        & !< net downward SW flux profile,
      &  flx_lw_net_clr(klev+1),    & !< clrsky downward LW flux profile,
      &  flx_sw_net_clr(klev+1),    & !< clrsky downward SW flux profile,
      &  flx_uplw_sfc(kbdim),             & !< sfc LW upward flux,
      &  flx_upsw_sfc(kbdim),             & !< sfc SW upward flux,
      &  flx_uplw_sfc_clr(kbdim),         & !< clrsky sfc LW upward flux,
      &  flx_upsw_sfc_clr(kbdim)            !< clrsky sfc SW upward flux,

    REAL(wp), INTENT(out), OPTIONAL ::    &
      &  flx_dnsw_diff_sfc(kbdim),        & !< sfc SW diffuse downward flux,
      &  flx_upsw_toa(kbdim),             & !< TOA SW upward flux,
      &  flx_dnpar_sfc(kbdim),            & !< PAR downward sfc flux
      &  vis_frc_sfc(kbdim),              & !< Visible fraction of net surface SW radiation
      &  nir_dff_frc_sfc(kbdim),          & !< Diffuse fraction of downward surface near-infrared radiation at surface
      &  vis_dff_frc_sfc(kbdim),          & !< Diffuse fraction of downward surface visible radiation at surface
      &  par_dff_frc_sfc(kbdim)             !< Diffuse fraction of downward surface PAR

    INTEGER  :: jk, jl, jp, jkb, jspec,   & !< loop indicies
      &  icldlyr(klev)                !< index for clear or cloudy

    REAL(wp) ::                           &
      &  zsemiss(jpband),           & !< LW surface emissivity by band
      &  ppd_hl(klev),              & !< pressure thickness in Pa
      &  pm_sfc(kbdim),                   & !< surface pressure in hPa
      &  amm,                             & !< molecular weight of moist air
      &  delta,                           & !< pressure thickness
      &  zscratch                           !< scratch array

    REAL(wp) :: z_sum_aea, z_sum_aes !help variables for aerosol

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

    REAL(wp) :: tune_dust(jpband)  ! local variable for LW absorption tuning of dust

    REAL(wp) ::                        &
         re_drop   (klev),       & !< effective radius of liquid
         re_cryst  (klev),       & !< effective radius of ice
         aux_out   (9),          &
         zmu0      (kbdim),            &
         zdayfrc   (kbdim)

    INTEGER, PARAMETER    :: rng_seed_size = 4
    INTEGER :: rnseeds(rng_seed_size)

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

    !
    ! 1.0 Constituent properties
    !--------------------------------

    !
    ! --- control for infintesimal cloud fractions
    !
    DO jk = 1, klev
      !
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

  tune_dust(1:jpband) = 1._wp
  DO jspec=1,jpband
    DO jk=1,klev
      jkb = klev+1-jk
      ! LW opt thickness of aerosols
      aer_tau_lw_vr(jk,jspec) =  zaeq1(jkb) * zaea_rrtm(jspec,1) &
      &                         + zaeq2(jl,jkb) * zaea_rrtm(jspec,2) &
      &   + tune_dust(jl,jspec) * zaeq3(jl,jkb) * zaea_rrtm(jspec,3) &
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

  CALL newcld_optics(                                                       &
  & klev         ,jpband       ,jpsw                                     ,&
  & zglac        ,zland        ,ktype        ,icldlyr      ,tk_fl_vr     ,&
  & zlwp_vr      ,ziwp_vr      ,zlwc_vr      ,ziwc_vr      ,cdnc_vr      ,&
  & cld_tau_lw_vr,cld_tau_sw_vr,cld_piz_sw_vr,cld_cg_sw_vr                )

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
  &    flx_dnsw        ,flx_upsw        ,flx_dnsw_clr    ,flx_upsw_clr,   &
  !    optional output
  &    flxd_dff_sfc=flx_dnsw_diff_sfc ,                                   &
  &    flxd_par_sfc    = flx_dnpar_sfc,                                   &
  &    vis_frc_sfc     = vis_frc_sfc,                                     &
  &    nir_dff_frc_sfc = nir_dff_frc_sfc,                                 &
  &    vis_dff_frc_sfc = vis_dff_frc_sfc,                                 &
  &    par_dff_frc_sfc = par_dff_frc_sfc                                  )


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
