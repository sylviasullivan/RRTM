  SUBROUTINE rrtm_interface(                                              &
    ! input
    & current_date                                                       ,&
    & jg              ,jb              ,irad                             ,&
    & jce             ,kbdim           ,klev                             ,&
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
    ! optional input
    & dust_tunefac                                                       ,&
    ! optional output
    & flx_dnsw_diff_sfc, flx_upsw_toa  ,flx_dnpar_sfc                    ,&
    & vis_frc_sfc     ,nir_dff_frc_sfc ,vis_dff_frc_sfc ,par_dff_frc_sfc  )

    TYPE(datetime), POINTER, INTENT(in) :: current_date
    
    INTEGER,INTENT(in)  ::                &
      &  jg,                              & !< domain index
      &  jb,                              & !< block index
      &  irad,                            & !< option for radiation scheme (RRTM/PSRAD); active in NWP mode only, 
      !                                        ECHAM mode uses a completely different interface
      &  jce,                             & !< number of columns
      &  kbdim,                           & !< first dimension of 2-d arrays
      &  klev                               !< number of levels

    INTEGER,INTENT(in)  ::                &
      &  ktype(kbdim)                       !< type of convection

    REAL(wp),INTENT(in) ::                &
      &  zland(kbdim),                    & !< land-sea mask. (1. = land, 0. = sea/lakes)
      &  zglac(kbdim),                    & !< fraction of land covered by glaciers
      &  pmu0(kbdim),                     & !< mu0 for solar zenith angle
      &  alb_vis_dir(kbdim),              & !< surface albedo for vis range and dir light
      &  alb_nir_dir(kbdim),              & !< surface albedo for NIR range and dir light
      &  alb_vis_dif(kbdim),              & !< surface albedo for vis range and dif light
      &  alb_nir_dif(kbdim),              & !< surface albedo for NIR range and dif light
      &  emis_rad(kbdim),                 & !< longwave surface emissivity
      &  pp_fl(kbdim,klev),               & !< full level pressure in Pa
      &  pp_hl(kbdim,klev+1),             & !< half level pressure in Pa
      &  pp_sfc(kbdim),                   & !< surface pressure in Pa
      &  tk_fl(kbdim,klev),               & !< full level temperature in K
      &  tk_hl(kbdim,klev+1),             & !< half level temperature in K
      &  tk_sfc(kbdim),                   & !< surface temperature in K
      &  xm_vap(kbdim,klev),              & !< specific humidity in g/g
      &  xm_liq(kbdim,klev),              & !< specific liquid water content
      &  xm_ice(kbdim,klev),              & !< specific ice content in g/g
      &  cdnc(kbdim,klev),                & !< cloud nuclei concentration
      &  cld_frc(kbdim,klev),             & !< fractional cloud cover
      &  xm_o3(kbdim,klev),               & !< o3 mass mixing ratio
      &  xm_co2(kbdim,klev),              & !< co2 mass mixing ratio
      &  xm_ch4(kbdim,klev),              & !< ch4 mass mixing ratio
      &  xm_n2o(kbdim,klev),              & !< n2o mass mixing ratio
      &  xm_cfc11(kbdim,klev),            & !< cfc 11 volume mixing ratio
      &  xm_cfc12(kbdim,klev),            & !< cfc 12 volume mixing ratio
      &  xm_o2(kbdim,klev),               & !< o2  mass mixing ratio
      &  zaeq1(kbdim,klev),               & !< aerosol continental
      &  zaeq2(kbdim,klev),               & !< aerosol maritime
      &  zaeq3(kbdim,klev),               & !< aerosol urban
      &  zaeq4(kbdim,klev),               & !< aerosol volcano ashes
      &  zaeq5(kbdim,klev)                  !< aerosol stratospheric background


    REAL(wp), INTENT(out) ::              &
      &  flx_lw_net(kbdim,klev+1),        & !< net downward LW flux profile,
      &  flx_sw_net(kbdim,klev+1),        & !< net downward SW flux profile,
      &  flx_lw_net_clr(kbdim,klev+1),    & !< clrsky downward LW flux profile,
      &  flx_sw_net_clr(kbdim,klev+1),    & !< clrsky downward SW flux profile,
      &  flx_uplw_sfc(kbdim),             & !< sfc LW upward flux,
      &  flx_upsw_sfc(kbdim),             & !< sfc SW upward flux,
      &  flx_uplw_sfc_clr(kbdim),         & !< clrsky sfc LW upward flux,
      &  flx_upsw_sfc_clr(kbdim)            !< clrsky sfc SW upward flux,

    REAL(wp), INTENT(in),  OPTIONAL ::    dust_tunefac(kbdim,jpband) ! LW absorption tuning factor for dust

    REAL(wp), INTENT(out), OPTIONAL ::    &
      &  flx_dnsw_diff_sfc(kbdim),        & !< sfc SW diffuse downward flux,
      &  flx_upsw_toa(kbdim),             & !< TOA SW upward flux,
      &  flx_dnpar_sfc(kbdim),            & !< PAR downward sfc flux
      &  vis_frc_sfc(kbdim),              & !< Visible fraction of net surface SW radiation
      &  nir_dff_frc_sfc(kbdim),          & !< Diffuse fraction of downward surface near-infrared radiation at surface
      &  vis_dff_frc_sfc(kbdim),          & !< Diffuse fraction of downward surface visible radiation at surface
      &  par_dff_frc_sfc(kbdim)             !< Diffuse fraction of downward surface PAR

    INTEGER  :: jk, jl, jp, jkb, jspec,   & !< loop indicies
      &  icldlyr(kbdim,klev)                !< index for clear or cloudy

    REAL(wp) ::                           &
      &  zsemiss(kbdim,jpband),           & !< LW surface emissivity by band
      &  ppd_hl(kbdim,klev),              & !< pressure thickness in Pa
      &  pm_sfc(kbdim),                   & !< surface pressure in hPa
      &  amm,                             & !< molecular weight of moist air
      &  delta,                           & !< pressure thickness
      &  zscratch                           !< scratch array

    REAL(wp) :: z_sum_aea, z_sum_aes !help variables for aerosol

    !
    ! --- vertically reversed _vr variables
    !
    REAL(wp) ::                           &
      &  col_dry_vr(kbdim,klev),          & !< number of molecules/cm2 of
      &  pm_fl_vr(kbdim,klev),            & !< full level pressure [hPa]
      &  pm_hl_vr(kbdim,klev+1),          & !< half level pressure [hPa]
      &  tk_fl_vr(kbdim,klev),            & !< full level temperature [K]
      &  tk_hl_vr(kbdim,klev+1),          & !< half level temperature [K]
      &  cdnc_vr(kbdim,klev),             & !< cloud nuclei concentration
      &  cld_frc_vr(kbdim,klev),          & !< secure cloud fraction
      &  ziwgkg_vr(kbdim,klev),           & !< specific ice water content
      &  ziwc_vr(kbdim,klev),             & !< ice water content per volume
      &  ziwp_vr(kbdim,klev),             & !< ice water path in g/m2
      &  zlwgkg_vr(kbdim,klev),           & !< specific liquid water content
      &  zlwp_vr(kbdim,klev),             & !< liquid water path in g/m2
      &  zlwc_vr(kbdim,klev),             & !< liquid water content per
      &  wkl_vr(kbdim,jpinpx,klev),       & !< number of molecules/cm2 of
      &  wx_vr(kbdim,jpxsec,klev),        & !< number of molecules/cm2 of
      &  cld_tau_lw_vr(kbdim,klev,jpband),& !< LW optical thickness of clouds
      &  cld_tau_sw_vr(kbdim,jpsw,klev),  & !< extincion
      &  cld_cg_sw_vr(kbdim,jpsw,klev),   & !< asymmetry factor
      &  cld_piz_sw_vr(kbdim,jpsw,klev),  & !< single scattering albedo
      &  aer_tau_lw_vr(kbdim,klev,jpband),& !< LW optical thickness of aerosols
      &  aer_tau_sw_vr(kbdim,klev,jpsw),  & !< aerosol optical thickness
      &  aer_cg_sw_vr(kbdim,klev,jpsw),   & !< aerosol asymmetry factor
      &  aer_piz_sw_vr(kbdim,klev,jpsw),  & !< aerosol single scattering albedo
      &  flx_uplw_vr(kbdim,klev+1),       & !< upward flux, total sky
      &  flx_uplw_clr_vr(kbdim,klev+1),   & !< upward flux, clear sky
      &  flx_dnlw_vr(kbdim,klev+1),       & !< downward flux, total sky
      &  flx_dnlw_clr_vr(kbdim,klev+1),   & !< downward flux, clear sky
      &  flx_upsw(kbdim,klev+1),          & !< upward flux total sky
      &  flx_upsw_clr(kbdim,klev+1),      & !< upward flux clear sky
      &  flx_dnsw(kbdim,klev+1),          & !< downward flux total sky
      &  flx_dnsw_clr(kbdim,klev+1)         !< downward flux clear sky

    REAL(wp) :: tune_dust(kbdim,jpband)  ! local variable for LW absorption tuning of dust

    ! Additional fields needed for PSRAD call
    LOGICAL ::                         &
         laland(kbdim),                & !< land sea mask, land=.true.
         laglac(kbdim)                   !< glacier mask, glacier=.true.

    REAL(wp) ::                        &
         re_drop   (kbdim,klev),       & !< effective radius of liquid
         re_cryst  (kbdim,klev),       & !< effective radius of ice
         aux_out   (kbdim,9),          &
         zmu0      (kbdim),            &
         zdayfrc   (kbdim)

    INTEGER, PARAMETER    :: rng_seed_size = 4
    INTEGER :: rnseeds(kbdim,rng_seed_size)

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

    IF (ltimer) CALL timer_start(timer_rrtm_prep)

    !
    ! --- control for infintesimal cloud fractions
    !
    DO jk = 1, klev
      !
      jkb = klev+1-jk
      cld_frc_vr(1:jce,jk)  = cld_frc(1:jce,jkb)

      DO jl=1,jce
        IF (cld_frc_vr(jl,jk) > 2.0_wp*EPSILON(1.0_wp)) THEN
          ! only clouds > 2 epsilon are made visible to radiation
          icldlyr  (jl,jk) = 1
          ziwgkg_vr(jl,jk) = xm_ice(jl,jkb)*1000.0_wp/cld_frc_vr(jl,jk)
          zlwgkg_vr(jl,jk) = xm_liq(jl,jkb)*1000.0_wp/cld_frc_vr(jl,jk)
        ELSE
          ! clouds <= 2 epsilon are ade invisble to radiation
          icldlyr  (jl,jk) = 0
          ziwgkg_vr(jl,jk) = 0.0_wp
          zlwgkg_vr(jl,jk) = 0.0_wp
        ENDIF
      END DO
    END DO
    !
    ! --- main constituent reordering
    !
    DO jl = 1, jce
      pm_hl_vr(jl,klev+1) = 0.01_wp*pp_hl(jl,1)
      tk_hl_vr(jl,klev+1) = tk_hl(jl,1)
      pm_sfc(jl)          = 0.01_wp*pp_sfc(jl)
    END DO

    DO jk = 1, klev
      jkb = klev+1-jk
      ! initialization
      wkl_vr(:,:,jk) = 0.0_wp
      wx_vr (:,:,jk) = 0.0_wp
      DO jl = 1, jce
        delta = pp_hl(jl,jkb+1)-pp_hl(jl,jkb)
        !
        ! --- thermodynamic arrays
        !
        pm_hl_vr(jl,jk) = 0.01_wp*pp_hl(jl,jkb+1)
        pm_fl_vr(jl,jk) = 0.01_wp*pp_fl(jl,jkb)
        tk_hl_vr(jl,jk) = tk_hl(jl,jkb+1)
        tk_fl_vr(jl,jk) = tk_fl(jl,jkb)
        !
        ! --- cloud properties
        !
        zscratch       = pp_fl(jl,jkb)/tk_fl(jl,jkb)
        ziwc_vr(jl,jk) = ziwgkg_vr(jl,jk)*zscratch/rd
        ziwp_vr(jl,jk) = ziwgkg_vr(jl,jk)*delta/grav
        zlwc_vr(jl,jk) = zlwgkg_vr(jl,jk)*zscratch/rd
        zlwp_vr(jl,jk) = zlwgkg_vr(jl,jk)*delta/grav
        cdnc_vr(jl,jk) = cdnc(jl,jkb)*1.e-6_wp
        !
        ! --- radiatively active gases
        !
        wkl_vr(jl,1,jk)   = xm_vap(jl,jkb)*amd/amw
        wkl_vr(jl,2,jk)   = xm_co2(jl,jkb)*amd/amco2
        wkl_vr(jl,3,jk)   = xm_o3(jl,jkb) *amd/amo3
        wkl_vr(jl,4,jk)   = xm_n2o(jl,jkb)*amd/amn2o
        wkl_vr(jl,6,jk)   = xm_ch4(jl,jkb)*amd/amch4
        wkl_vr(jl,7,jk)   = xm_o2 (jl,jkb)*amd/amo2
        amm               = (1.0_wp-wkl_vr(jl,1,jk))*amd + wkl_vr(jl,1,jk)*amw
        col_dry_vr(jl,jk) = (0.01_wp*delta)*10.0_wp*avo/grav/amm / (1.0_wp+wkl_vr(jl,1,jk))
        !
        ! --- alternate treatment for cfcs
        !
        wx_vr(jl,2,jk) = col_dry_vr(jl,jk)*xm_cfc11(jl,jkb)*1.e-20_wp
        wx_vr(jl,3,jk) = col_dry_vr(jl,jk)*xm_cfc12(jl,jkb)*1.e-20_wp
      END DO
    END DO
    DO jp = 1, 7
      wkl_vr(1:jce,jp,:)=col_dry_vr(1:jce,:)*wkl_vr(1:jce,jp,:)
    END DO
    !
    ! 2.0 Surface Properties
    ! --------------------------------
    zsemiss(1:jce,:) = SPREAD(emis_rad(1:jce),2,jpband)
    !
    ! 3.0 Particulate Optical Properties
    ! --------------------------------
    ppd_hl(1:jce,:) = pp_hl(1:jce,2:klev+1)-pp_hl(1:jce,1:klev)

    IF (PRESENT(dust_tunefac)) THEN
       tune_dust(1:jce,1:jpband) = dust_tunefac(1:jce,1:jpband)
    ELSE
      tune_dust(1:jce,1:jpband) = 1._wp
    ENDIF
    DO jspec=1,jpband
      DO jk=1,klev
        jkb = klev+1-jk
        DO jl = 1,jce
          ! LW opt thickness of aerosols
          aer_tau_lw_vr(jl,jk,jspec) =  zaeq1(jl,jkb) * zaea_rrtm(jspec,1) &
            &                         + zaeq2(jl,jkb) * zaea_rrtm(jspec,2) &
            &   + tune_dust(jl,jspec) * zaeq3(jl,jkb) * zaea_rrtm(jspec,3) &
            &                         + zaeq4(jl,jkb) * zaea_rrtm(jspec,4) &
            &                         + zaeq5(jl,jkb) * zaea_rrtm(jspec,5)
        ENDDO
      ENDDO
    ENDDO
    DO jspec=1+jpband,jpband+jpsw
      DO jk=1,klev
        jkb = klev+1-jk
        DO jl = 1,jce

          z_sum_aea = zaeq1(jl,jkb) * zaea_rrtm(jspec,1) &
            &       + zaeq2(jl,jkb) * zaea_rrtm(jspec,2) &
            &       + zaeq3(jl,jkb) * zaea_rrtm(jspec,3) &
            &       + zaeq4(jl,jkb) * zaea_rrtm(jspec,4) &
            &       + zaeq5(jl,jkb) * zaea_rrtm(jspec,5)

          z_sum_aes = zaeq1(jl,jkb) * zaes_rrtm(jspec,1) &
            &       + zaeq2(jl,jkb) * zaes_rrtm(jspec,2) &
            &       + zaeq3(jl,jkb) * zaes_rrtm(jspec,3) &
            &       + zaeq4(jl,jkb) * zaes_rrtm(jspec,4) &
            &       + zaeq5(jl,jkb) * zaes_rrtm(jspec,5)

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
    ENDDO
    
    IF (lrad_aero_diag) THEN
      CALL rad_aero_diag (                                  &
      & jg              ,jb              ,jce             , &
      & kbdim           ,klev            ,jpband          , &
      & jpsw            ,aer_tau_lw_vr   ,aer_tau_sw_vr   , &
      & aer_piz_sw_vr   ,aer_cg_sw_vr                       )
    END IF

    IF (irad == 1) THEN
      CALL newcld_optics(                                                       &
        & jce          ,kbdim        ,klev         ,jpband       ,jpsw         ,&
        & zglac        ,zland        ,ktype        ,icldlyr      ,tk_fl_vr     ,&
        & zlwp_vr      ,ziwp_vr      ,zlwc_vr      ,ziwc_vr      ,cdnc_vr      ,&
        & cld_tau_lw_vr,cld_tau_sw_vr,cld_piz_sw_vr,cld_cg_sw_vr                )
    ELSE
      DO jl = 1,jce
        IF (zland(jl) >= 0.5_wp) THEN
          laland(jl) = .TRUE.
        ELSE
          laland(jl) = .FALSE.
        ENDIF
        IF (zglac(jl) >= 0.5_wp) THEN
          laglac(jl) = .TRUE.
        ELSE
          laglac(jl) = .FALSE.
        ENDIF
      ENDDO
      CALL psrad_cloud_optics(                                          &
         & laglac        ,laland        ,jce           ,kbdim          ,& 
         & klev          , ktype        ,&
         & icldlyr       ,zlwp_vr       ,ziwp_vr       ,zlwc_vr        ,&
         & ziwc_vr       ,cdnc_vr       ,cld_tau_lw_vr ,cld_tau_sw_vr  ,&
         & cld_piz_sw_vr ,cld_cg_sw_vr  ,re_drop       ,re_cryst    )  
    ENDIF

    IF (ltimer) CALL timer_stop(timer_rrtm_prep)

    !
    ! 4.0 Radiative Transfer Routines
    ! --------------------------------
    IF (ltimer) CALL timer_start(timer_lrtm)
    IF (irad == 1) THEN
      CALL lrtm(                                                                &
        !    input
        &    jce             ,klev                                             ,&
        &    pm_fl_vr        ,pm_sfc          ,tk_fl_vr        ,tk_hl_vr       ,&
        &    tk_sfc          ,wkl_vr          ,wx_vr           ,col_dry_vr     ,&
        &    zsemiss         ,cld_frc_vr      ,cld_tau_lw_vr   ,aer_tau_lw_vr  ,&
        !    output
        &    flx_uplw_vr     ,flx_dnlw_vr     ,flx_uplw_clr_vr,flx_dnlw_clr_vr )
    ELSE
      ! Seeds for random numbers come from least significant digits of pressure field 
      !
      rnseeds(1:jce,1:rng_seed_size) = (pm_fl_vr(1:jce,1:rng_seed_size) -  &
         int(pm_fl_vr(1:jce,1:rng_seed_size)))* 1E9
      CALL psrad_lrtm(jce                                                       ,&
           & kbdim           ,klev            ,pm_fl_vr        ,pm_sfc          ,&
           & tk_fl_vr        ,tk_hl_vr        ,tk_sfc          ,wkl_vr          ,&
           & wx_vr           ,col_dry_vr      ,zsemiss         ,cld_frc_vr      ,&
           & cld_tau_lw_vr   ,aer_tau_lw_vr   ,rnseeds         ,flx_uplw_vr     ,&     
           & flx_dnlw_vr     ,flx_uplw_clr_vr ,flx_dnlw_clr_vr )
    ENDIF
    IF (ltimer) CALL timer_stop(timer_lrtm)


    IF (ltimer) CALL timer_start(timer_srtm)
    IF (irad == 1) THEN
      CALL srtm_srtm_224gp(                                                     &
        !    input
        &    jce             ,kbdim           ,klev            ,jpsw           ,&
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
    ELSE
      ! Reset random seeds so SW doesn't depend on what's happened in LW but is also independent
      !
      rnseeds(1:jce,1:rng_seed_size) = (pm_fl_vr(1:jce,rng_seed_size:1:-1) - &
         int(pm_fl_vr(1:jce,rng_seed_size:1:-1)))* 1E9
      WHERE (pmu0(1:jce) > 0.0_wp)
         zdayfrc(1:jce) = 1.0_wp
      ELSEWHERE
         zdayfrc(1:jce) = 0.0_wp
      END WHERE
      zmu0(1:jce) = MAX(pmu0(1:jce),0.05_wp)
      
      !
      CALL psrad_srtm(jce                                                      , & 
         &  kbdim           ,klev            ,pm_fl_vr        ,tk_fl_vr        , &
         &  wkl_vr          ,col_dry_vr      ,alb_vis_dir     ,alb_vis_dif     , &
         &  alb_nir_dir     ,alb_nir_dif     ,zmu0, zdayfrc   ,ssi_radt/psctm  , &
         &  psctm           ,cld_frc_vr      ,cld_tau_sw_vr   ,cld_cg_sw_vr    , &
         &  cld_piz_sw_vr   ,aer_tau_sw_vr   ,aer_cg_sw_vr    ,aer_piz_sw_vr   , & 
         &  rnseeds         ,flx_dnsw        ,flx_upsw        ,flx_dnsw_clr    , &
         & flx_upsw_clr     ,aux_out(:,1)    ,aux_out(:,2)    ,aux_out(:,3)    , &
         &  aux_out(:,4)    ,aux_out(:,5)    ,aux_out(:,6)    ,aux_out(:,7)    , &
         aux_out(:,8)    ,aux_out(:,9) )

      !   dnpar_sfc        = dnpar_sfc_dir    + dnpar_sfc_dif                 
      flx_dnpar_sfc(1:jce) = aux_out(1:jce,2) + aux_out(1:jce,5)

      ! Reset solar fluxes to zero at dark points
      DO jl = 1, jce
        IF (pmu0(jl) <= 0._wp) THEN
          flx_dnsw(jl,:) = 0._wp
          flx_upsw(jl,:) = 0._wp
        ENDIF
      ENDDO
    ENDIF
    IF (ltimer) CALL timer_stop(timer_srtm)


    ! 5.0 Post Processing
    ! --------------------------------
    IF (ltimer) CALL timer_start(timer_rrtm_post)

    DO jk = 1, klev+1
      jkb = klev+2-jk
      DO jl = 1, jce
        flx_lw_net(jl,jk)     = flx_dnlw_vr(jl,jkb)-flx_uplw_vr(jl,jkb)
        flx_lw_net_clr(jl,jk) = flx_dnlw_clr_vr(jl,jkb)-flx_uplw_clr_vr(jl,jkb)
        flx_sw_net(jl,jk)     = flx_dnsw(jl,jk) - flx_upsw(jl,jk)
        flx_sw_net_clr(jl,jk) = flx_dnsw_clr(jl,jk)-flx_upsw_clr(jl,jk)
      END DO
    END DO
    flx_uplw_sfc(1:jce)     = flx_uplw_vr(1:jce,1)
    flx_uplw_sfc_clr(1:jce) = flx_uplw_clr_vr(1:jce,1)
    flx_upsw_sfc(1:jce)     = flx_upsw(1:jce,klev+1)
    flx_upsw_sfc_clr(1:jce) = flx_upsw_clr(1:jce,klev+1)
    IF (PRESENT(flx_upsw_toa)) flx_upsw_toa(1:jce) = flx_upsw(1:jce,1)
    IF (irad /= 1 .AND. PRESENT(flx_dnsw_diff_sfc))    &  ! approximate calculation!!
      !   dnsw_diff_sfc        = vis_dn_dff_sfc   + nir_dn_dff_sfc
      flx_dnsw_diff_sfc(1:jce) = aux_out(1:jce,4) + aux_out(1:jce,6)
!!$    sw_irr_toa(1:jce)       = flx_dnsw(1:jce,1)
    !
    IF (ltimer) CALL timer_stop(timer_rrtm_post)

  END SUBROUTINE rrtm_interface
