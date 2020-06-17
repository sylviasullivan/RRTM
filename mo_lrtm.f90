  !-----------------------------------------------------------------------------
  !>
  !! @brief Prepares information for radiation call
  !!
  !! @remarks: This program is the driver subroutine for RRTMG_LW, the AER LW
  !! radiation model for application to GCMs, that has been adapted from RRTM_LW
  !! for improved efficiency.
  !! This routine:
  !!    1) calls INATM to read in the atmospheric profile from GCM;
  !!       all layering in RRTMG is ordered from surface to toa.
  !! >> 20200306, sylvia The contents of this subroutine have been pasted here.
  !!    2) calls COEFFS to calculate various quantities needed for
  !!       the radiative transfer algorithm
  !!    3) calls TAUMOL to calculate gaseous optical depths for each
  !!       of the 16 spectral bands
  !!    4) calls RTRNMR (for both clear and cloudy profiles) to perform the
  !!       radiative transfer calculation with a maximum-random cloud
  !!       overlap method, or calls RTRN to use random cloud overlap.
  !!    5) passes the necessary fluxes and cooling rates back to GCM
  !!
  !! The optional calculation of the change in upward flux as a function of
  !! surface temperature is available (controlled by input flag idrv).  This
  !! can be utilized  to approximate adjustments to the upward flux profile
  !! caused only by a change in surface temperature between full radiation
  !! calls.  This feature uses the pre-calculated derivative of the Planck
  !! function with respect to surface temperature. To calculate these
  !! derivatives idrv must be set to 1.
  !'
! >> sylvia_20200318
! Return to making this a module not a subroutine.
MODULE mo_lrtm

  ! >> sylvia_20200318
  ! Link to other modules.
  USE mo_lrtm_taumol, ONLY: lrtm_taumol
  USE mo_lrtm_coeffs, ONLY: lrtm_coeffs
  USE mo_lrtm_rtrnmr, ONLY: lrtm_rtrnmr

CONTAINS

  SUBROUTINE lrtm     (                            &
       & klev       ,play     ,psfc               ,&
       & tlay       ,tlev     ,tsfc     ,wkl_2d   ,&
       & wx_2d      ,coldry_2d,emis     ,cldfr    ,&
       & taucld     ,tauaer   ,uflx     ,dflx     ,&
       & uflxc      ,dflxc    ,duflx_dt ,duflxc_dt)

    ! >> sylvia_20200310
    ! Set up the wp type that is often used below.
    ! Pull the nbndlw and ngptlw constants from mo_lrtm_par.f90
    ! >> sylvia_20200313
    ! Add amwm, amd variables from shared/mo_physical_constants.f90
    ! Add grav from mo_physical_constants
    ! >> sylvia_20200318
    ! Add nmol from mo_lrtm_par.f90
    INTEGER, PARAMETER :: pd =  12
    INTEGER, PARAMETER :: rd = 307
    INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(pd,rd) !< double precision
    INTEGER, PARAMETER :: wp = dp
    INTEGER, PARAMETER :: nbndlw = 16
    INTEGER, PARAMETER :: ngptlw = 140 ! set to 256 for 256 gpt model
    REAL(wp), PARAMETER :: amw = 18.0154_wp       !! [g/mol] H2O
    REAL(wp), PARAMETER :: amd = 28.970_wp        !> [g/mol] dry air
    REAL(wp), PARAMETER :: grav = 9.80665_wp     !> [m/s2] av. gravitational acceleration
    INTEGER, PARAMETER :: nmol = 7
    ! << sylvia_20200310
    
    INTEGER, INTENT(in) :: &
         & klev             !< Number of model layers

    REAL(wp), INTENT(in) :: &
         & play(:),       & !< Layer pressures [hPa, mb] (klev)
         & psfc,          & !< Surface pressure [hPa, mb]
         & tlay(:),       & !< Layer temperatures [K] (klev)
         & tlev(:),       & !< Interface temperatures [K] (klev+1)
         & tsfc,          & !< Surface temperature [K] 
         & wkl_2d(:,:),   & !< Gas volume mixing ratios
         & wx_2d(:,:),    & !< CFC type gas volume mixing ratios
         & coldry_2d(:),  & !< Column dry amount
         & emis(:),       & !< Surface emissivity  (nbndlw)
         & cldfr(:),      & !< Cloud fraction  (klev)
         & taucld(:,:),   & !< Coud optical depth (klev,nbndlw)
         & tauaer(:,:)      !< Qerosol optical depth (klev,nbndlw)

    REAL(wp), INTENT(out) :: &
         uflx(:) ,  & !< Tot sky longwave upward flux [W/m2], (klev+1)
         dflx(:) ,  & !< Tot sky longwave downward flux [W/m2], (klev+1)
         uflxc(:),  & !< Clr sky longwave upward flux [W/m2], (klev+1)
         dflxc(:)     !< Clr sky longwave downward flux [W/m2], (klev+1)

    REAL(wp), INTENT(out), OPTIONAL :: &
         & duflx_dt(:),  & !< change in uflx wrt Temp [w/m2/k] (klev)
         & duflxc_dt(:)    !< change in uflxc wrt Temp [w/m2/k] (klev)

    INTEGER :: ncbands     !< number of cloud spectral bands

    REAL(wp) ::               &
         & tz(0:klev),           & !< level (interface) temperatures [K]
         & wbrodl(klev),         & !< broadening gas column density (mol/cm2)
         & pwvcm,                & !< precipitable water vapor [cm]
         & fracs(klev,ngptlw),   & !< layer cloud fraction
         & taug(klev,ngptlw),    & !< gas optical depth
         & taut(klev,ngptlw),    & !< gaseous + aerosol optical depths
         & taucloud(klev,nbndlw),& !< layer in-cloud optical depth
         & totuflux(0:klev),     & !< upward longwave flux (w/m2)
         & totdflux(0:klev),     & !< downward longwave flux (w/m2)
         & fnet(0:klev),         & !< net longwave flux (w/m2)
         & totuclfl(0:klev),     & !< clear sky upward longwave flux (w/m2)
         & totdclfl(0:klev),     & !< clear sky downward longwave flux (w/m2)
         & fnetc(0:klev),        & !< clear sky net longwave flux (w/m2)
         & dtotuflux_dt(0:klev), & !< change in uflxx (w/m2/k) wrt Temp
         & dtotuclfl_dt(0:klev), & !< change in cuflx  (w/m2/k) wrt Temp
         & tauctot(klev)             !< band integrated cloud optical depth

    INTEGER ::            &
         & laytrop,              & !< tropopause layer index
         & jp(klev),             & !< lookup table index
         & jt(klev),             & !< lookup table index
         & jt1(klev),            & !< lookup table index
         & indself(klev),        &
         & indfor(klev),         &
         & indminor(klev)

    REAL(wp) ::                  &
         & planklay(klev,nbndlw),   & !
         & planklev(0:klev,nbndlw), & !
         & plankbnd(nbndlw),        & !
         & dplankbnd_dt(nbndlw),    & !
         & colh2o(klev),            & !< column amount (h2o)
         & colco2(klev),            & !< column amount (co2)
         & colo3(klev),             & !< column amount (o3)
         & coln2o(klev),            & !< column amount (n2o)
         & colco(klev),             & !< column amount (co)
         & colch4(klev),            & !< column amount (ch4)
         & colo2(klev),             & !< column amount (o2)
         & colbrd(klev),            & !< column amount (broadening gases)
         & selffac(klev),           &
         & selffrac(klev),          &
         & forfac(klev),            &
         & forfrac(klev),           &
         & minorfrac(klev),         &
         & scaleminor(klev),        &
         & scaleminorn2(klev)

    REAL(wp) ::      &
         fac00(klev),        &
         fac01(klev),        &
         fac10(klev),        &
         fac11(klev),        &
         rat_h2oco2(klev),   &
         rat_h2oco2_1(klev), &
         rat_h2oo3(klev),    &
         rat_h2oo3_1(klev),  &
         rat_h2on2o(klev),   &
         rat_h2on2o_1(klev), &
         rat_h2och4(klev),   &
         rat_h2och4_1(klev), &
         rat_n2oco2(klev),   &
         rat_n2oco2_1(klev), &
         rat_o3co2(klev),    &
         rat_o3co2_1(klev)

    INTEGER :: jk, ib, ig ! loop indicies

    INTEGER,  PARAMETER :: idrv   =  0 !< Flag for calculating derivs
    INTEGER,  PARAMETER :: iout   =  0 !< option output flag
    INTEGER,  PARAMETER :: iend   = 16 !< last band
    INTEGER,  PARAMETER :: istart =  1 !< first band
    REAL(wp), PARAMETER :: cldmin = 1.e-20_wp ! minimum val for clouds

    ! >> sylvia_20200306
    ! Adding in variables for the inatm call, where they have INTENT(out).
    ! >> sylvia_20200310
    ! Removing tz(0:klev), wbrodl(:), pwvcm(:) below as they have already been declared
    ! >> sylvia_20200318
    ! Adding in ngb from mo_lrtm_par module.
    REAL(wp) :: amttl, wvttl, wvsh, summol
    INTEGER  :: imol         !< Loop indices
    INTEGER :: ngb(ngptlw)
    ! << sylvia_20200306
   
    ! >> sylvia_20200318
    ! Copying and pasting ngb values from mo_lrtm_setup.
    ngb(:) = (/1,1,1,1,1,1,1,1,1,1, &                 ! band 1
         2,2,2,2,2,2,2,2,2,2,2,2, &                   ! band 2
         3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3, &           ! band 3
         4,4,4,4,4,4,4,4,4,4,4,4,4,4, &               ! band 4
         5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5, &           ! band 5
         6,6,6,6,6,6,6,6, &                           ! band 6
         7,7,7,7,7,7,7,7,7,7,7,7, &                   ! band 7
         8,8,8,8,8,8,8,8, &                           ! band 8
         9,9,9,9,9,9,9,9,9,9,9,9, &                   ! band 9
         10,10,10,10,10,10, &                         ! band 10
         11,11,11,11,11,11,11,11, &                   ! band 11
         12,12,12,12,12,12,12,12, &                   ! band 12
         13,13,13,13, &                               ! band 13
         14,14, &                                     ! band 14
         15,15, &                                     ! band 15
         16,16/)
    ! << sylvia_20200318 
    !
    ! 1.0 Convert 2D data into 1D column data an prepare some auxilliary info
    ! --------------------------------
    !
    ! >> sylvia_20200306
    ! Copy-and-pasting the contents of subroutine inatm here.
    !  --- Initialize all molecular amounts and cloud properties to zero
    !
    amttl  = 0.0_wp
    wvttl  = 0.0_wp
    !
    !  --- Remap data to 1D column arrays
    !
    tz(0:klev) = tlev(1:klev+1)
    !
    !  --- Calculate auxillary information
    !
    DO jk = 1, klev
      summol = 0.0_wp
!CDIR EXPAND=nmol
      DO imol = 2, nmol
        summol = summol + wkl_2d(imol,jk)
      ENDDO
      wbrodl(jk) = coldry_2d(jk) - summol
      amttl = amttl + coldry_2d(jk)+wkl_2d(1,jk)
      wvttl = wvttl + wkl_2d(1,jk)
    ENDDO

    wvsh = (amw * wvttl) / (amd * amttl)
    pwvcm = wvsh * (1.e3_wp * psfc) / (1.e2_wp * grav)
    ! << sylvia_20200306
    
    ncbands     = 1
    tauctot(:)  = 0.0_wp

    DO ib = 1,nbndlw
      DO jk = 1, klev
        tauctot(jk)  = tauctot(jk) + taucld(jk,ib)
      END DO
    ENDDO

    DO jk = 1, klev
      IF (cldfr(jk) .GE. cldmin .AND. tauctot(jk) .GE. cldmin) THEN
        ncbands = 16
        ! >> sylvia_20200318
        ! Removing CDIR EXPAND=nbndlw from below.
        DO ib = 1,nbndlw
          taucloud(jk,ib) = taucld(jk,ib)
        END DO
      ELSE
        taucloud(jk,:) = 0.0_wp
      END IF
    ENDDO

    !
    ! 2.0  Calculate information needed by the radiative transfer routine
    ! that is specific to this atmosphere, especially some of the
    ! coefficients and indices needed to compute the optical depths
    ! by interpolating data from stored reference atmospheres.
    ! --------------------------------
    !
    CALL lrtm_coeffs(                                                &
         & klev         ,istart       ,play         ,tlay           ,&
         & tz           ,tsfc         ,emis         ,coldry_2d      ,&
         & wkl_2d       ,wbrodl       ,laytrop      ,jp             ,&
         & jt           ,jt1          ,planklay     ,planklev       ,&
         & plankbnd     ,idrv         ,dplankbnd_dt ,colh2o         ,&
         & colco2       ,colo3        ,coln2o       ,colco          ,&
         & colch4       ,colo2        ,colbrd       ,fac00          ,&
         & fac01        ,fac10        ,fac11        ,rat_h2oco2     ,&
         & rat_h2oco2_1 ,rat_h2oo3    ,rat_h2oo3_1  ,rat_h2on2o     ,&
         & rat_h2on2o_1 ,rat_h2och4   ,rat_h2och4_1 ,rat_n2oco2     ,&
         & rat_n2oco2_1 ,rat_o3co2    ,rat_o3co2_1  ,selffac        ,&
         & selffrac     ,indself      ,forfac       ,forfrac        ,&
         & indfor       ,minorfrac    ,scaleminor   ,scaleminorn2   ,&
         & indminor     )

    !
    !  3.0 Calculate the gaseous optical depths and Planck fractions for
    !  each longwave spectral band.
    ! --------------------------------
    !
    CALL lrtm_taumol(                                                &
         & klev         ,play         ,wx_2d        ,coldry_2d      ,&
         & laytrop      ,jp           ,jt           ,jt1            ,&
         & colh2o       ,colco2       ,colo3        ,coln2o         ,&
         & colco        ,colch4       ,colo2        ,colbrd         ,&
         & fac00        ,fac01        ,fac10        ,fac11          ,&
         & rat_h2oco2   ,rat_h2oco2_1 ,rat_h2oo3    ,rat_h2oo3_1    ,&
         & rat_h2on2o   ,rat_h2on2o_1 ,rat_h2och4   ,rat_h2och4_1   ,&
         & rat_n2oco2   ,rat_n2oco2_1 ,rat_o3co2    ,rat_o3co2_1    ,&
         & selffac      ,selffrac     ,indself      ,forfac         ,&
         & forfrac      ,indfor       ,minorfrac    ,scaleminor     ,&
         & scaleminorn2 ,indminor     ,fracs        ,taug)
    ! --- Combine gaseous and aerosol optical depths, if aerosol active
    !
    DO ig = 1, ngptlw
      DO jk = 1, klev
          taut(jk,ig) = taug(jk,ig) + tauaer(jk,ngb(ig))
      ENDDO
    ENDDO
    !
    ! 4.0 Call the radiative transfer routine. (random-maximum overlap)
    ! --------------------------------
    !
    CALL lrtm_rtrnmr(                                                &
         & klev         ,istart       ,iend         ,iout           ,&
         & emis         ,ncbands      ,cldfr        ,taucloud       ,&
         & planklay     ,planklev     ,plankbnd     ,pwvcm          ,&
         & fracs        ,taut         ,totuflux     ,totdflux       ,&
         & fnet         ,totuclfl     ,totdclfl     ,fnetc          ,&
         & idrv         ,dplankbnd_dt ,dtotuflux_dt ,dtotuclfl_dt   )
    !
    ! 5.0 Finalize output
    ! --------------------------------
    !
    DO jk = 0, klev
        uflx(jk+1)  = totuflux(jk)
        dflx(jk+1)  = totdflux(jk)
        uflxc(jk+1) = totuclfl(jk)
        dflxc(jk+1) = totdclfl(jk)
    ENDDO

    IF (idrv .EQ. 1 .and. present(duflx_dt)) THEN
      DO jk = 0, klev
          duflx_dt(jk+1) = dtotuflux_dt(jk)
      ENDDO
    ENDIF
    IF (idrv .EQ. 1 .and. present(duflxc_dt) ) THEN
      DO jk = 0, klev
          duflxc_dt(jk+1) = dtotuclfl_dt(jk)
      ENDDO
    END IF

  END SUBROUTINE lrtm

END MODULE mo_lrtm
