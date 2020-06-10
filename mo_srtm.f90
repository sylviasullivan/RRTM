! >> sylvia_20200319
! Combine everything into a module.
MODULE mo_srtm

    USE mo_srtm_config, ONLY : delwave, wavenum2, wavenum1,            &
    &    ngc, jpinpx, jpb1, jpb2, preflog, tref, repclc, replog,     &
    &    ssi_default
    USE mo_srtm_taumol, ONLY :                                         &
    &    srtm_taumol16, srtm_taumol17, srtm_taumol18, srtm_taumol19, &
    &    srtm_taumol20, srtm_taumol21, srtm_taumol22, srtm_taumol23, &
    &    srtm_taumol24, srtm_taumol25, srtm_taumol26, srtm_taumol27, &
    &    srtm_taumol28, srtm_taumol29

    ! >> sylvia_20200311
    ! Set up the wp type that is often used below.
    ! Pull the nbndlw and ngptlw constants from mo_lrtm_par.f90
    INTEGER, PARAMETER :: pd = 12
    INTEGER, PARAMETER :: rdd = 307
    INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(pd,rdd) !< double precision
    INTEGER, PARAMETER :: wp = dp
    INTEGER, PARAMETER :: pi4 =  9
    INTEGER, PARAMETER :: i4 = SELECTED_INT_KIND(pi4)   !< at least 4 byte integer
    ! << sylvia_20200311

    ! >> sylvia_20200311
    ! Pull the constants from mo_physical_constants.
    REAL(wp), PARAMETER :: rd    = 287.04_wp        !> [J/K/kg] gas constant
    REAL(wp), PARAMETER :: rv    = 461.51_wp        !> [J/K/kg] gas constant for water vapor
    REAL(wp), PARAMETER :: rdv   = rd/rv            !> [ ]
    REAL(wp), PARAMETER :: grav  = 9.80665_wp       !> [m/s2] av. gravitational acceleration
    REAL(wp), PARAMETER :: rgrav = 1._wp/grav       !! [s2/m]
    ! << sylvia_20200311

CONTAINS

  SUBROUTINE srtm_srtm_224gp                                                   &
    !  input
    & (klev            , ksw                                                 , &
    &  alb_vis_dir     , alb_nir_dir     , alb_vis_dif     , alb_nir_dif     , &
    &  pm_fl_vr        , tk_fl_vr        , prmu0                             , &
    &  col_dry_vr      , wkl_vr                                              , &
    &  cld_frc_vr      , cld_tau_sw_vr   , cld_cg_sw_vr    , cld_piz_sw_vr   , &
    &  aer_tau_sw_vr   , aer_cg_sw_vr    , aer_piz_sw_vr                     , &
    &  ssi                                                                   , &
    !  output
    &  flxd_sw         , flxu_sw         , flxd_sw_clr     , flxu_sw_clr     , &
    ! optional output
    &  flxd_dff_sfc    , flxd_par_sfc    , vis_frc_sfc                       , &
    &  nir_dff_frc_sfc , vis_dff_frc_sfc , par_dff_frc_sfc                     )
    
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!! PARAMETER BLOCK !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    ! >> sylvia_20200311
    ! jpinpx value taken from mo_srtm_config.
    INTEGER, PARAMETER   :: jpinpx = 35
    ! << sylvia_20200311
    
    REAL(wp), PARAMETER :: nir_vis_boundary   = 14500._wp
    REAL(wp), PARAMETER :: par_lower_boundary = 14285.7143_wp ! equivalent to
    !                                                         ! 700nm wavelength
    REAL(wp), PARAMETER :: par_upper_boundary = 25000._wp     ! equivalent to
    !                                                         ! 400nm wavelength

    ! >> sylvia_20200311
    ! Pulling these values from mo_srtm_config.f90
    REAL(wp), PARAMETER :: ssi_default(14) =  (/ & !< SRTM default solar flux (W/m2) in 14 SW bands
      & 12.1095682699999987_wp   , 20.3650825429849398_wp   , 23.7297328613475429_wp   , &
      & 22.4276934179066849_wp   , 55.6266126199999960_wp   , 1.0293153082385436E+02_wp, &
      & 24.2936128100154392_wp   , 3.4574251380000004E+02_wp, 2.1818712729860400E+02_wp, &
      & 3.4719231470000005E+02_wp, 1.2949501812000000E+02_wp, 50.1522503011060650_wp   , &
      & 3.0799387010047838_wp    , 12.8893773299999985_wp  /)
    ! sum of 14 bands is: 1.3682223735968237E+03
    ! << sylvia_20200311
    
    REAL(wp), PARAMETER :: zdecorr = 2000.0_wp ! decorrelation length scale del(z0)
  
    !-- Interface to RRTM_SW
    !     JJMorcrette 030225
    !     JJMorcrette 20071015 3D fields of CO2, CH4, N2O and NO2
    !        D.Salmond  31-Oct-2007 Vector version in the style of RRTM from Meteo France & NEC
    !     define bands 24-28 as visible
    !     TJ Raddatz  20091223 diagnosis of photosynthetically active radiation (PAR)
    !     TJ Raddatz  20100111 passing visible and NIR surface albedo to
    !        radiative transfer calculations
    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!! ENDOF PARAMETER BLOCK !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !-- Input arguments
    
    ! >> sylvia_20200311
    ! Removing kbdim as a dimension below.
    ! No INTENT(in) for alb_*_* type variables.
    INTEGER,  INTENT(in)    :: klev
    INTEGER,  INTENT(in)    :: ksw
    REAL(wp)                :: alb_vis_dir
    REAL(wp)                :: alb_nir_dir ! >> sylvia_20200310, intent(in) removed
    REAL(wp)                :: alb_vis_dif
    REAL(wp)                :: alb_nir_dif
    REAL(wp), INTENT(in)    :: prmu0
    REAL(wp), INTENT(in)    :: pm_fl_vr(klev)
    REAL(wp), INTENT(in)    :: tk_fl_vr(klev)
    REAL(wp), INTENT(in)    :: col_dry_vr(klev)
    REAL(wp), INTENT(in)    :: wkl_vr(jpinpx,klev)
    REAL(wp), INTENT(in)    :: cld_frc_vr(klev)
    REAL(wp), INTENT(in)    :: cld_tau_sw_vr(ksw,klev)
    REAL(wp), INTENT(in)    :: cld_cg_sw_vr(ksw,klev)
    REAL(wp), INTENT(in)    :: cld_piz_sw_vr(ksw,klev)
    ! hs: The order of indices in the aerosol optical properties was wrong. However,
    ! it should be checked if a fix in the interface would be more appropriate.
    !   REAL(wp), INTENT(in)    :: aer_tau_sw_vr(ksw,klev)
    !   REAL(wp), INTENT(in)    :: aer_cg_sw_vr(ksw,klev)
    !   REAL(wp), INTENT(in)    :: aer_piz_sw_vr(ksw,klev)
    REAL(wp), INTENT(in)    :: aer_tau_sw_vr(klev,ksw)
    REAL(wp), INTENT(in)    :: aer_cg_sw_vr(klev,ksw)
    REAL(wp), INTENT(in)    :: aer_piz_sw_vr(klev,ksw)
    REAL(wp), INTENT(in)    :: ssi(ksw)

    !-- Output arguments

    REAL(wp), INTENT(out)   :: flxd_sw(klev+1)     !< downward flux total sky
    REAL(wp), INTENT(out)   :: flxd_sw_clr(klev+1) !< downward flux clear sky
    REAL(wp), INTENT(out)   :: flxu_sw(klev+1)     !< upward flux total sky
    REAL(wp), INTENT(out)   :: flxu_sw_clr(klev+1) !< upward flux clear sky

    REAL(wp), INTENT(out), OPTIONAL :: &
      & flxd_dff_sfc,     & !< surface downward diffuse rad
      & flxd_par_sfc,     & !< surface downward photosynthetically active rad
      & vis_frc_sfc,      & !< Visible fraction of net surface radiation
      & nir_dff_frc_sfc,  & !< Diffuse fraction of downward surface near-infrared radiation
      & vis_dff_frc_sfc,  & !< Diffuse fraction of downward surface visible radiation
      & par_dff_frc_sfc     !< Diffuse fraction of downward surface PAR

    !-- Local variables

    REAL(wp) :: z_colmol(klev) ,  z_co2mult(klev)
    REAL(wp) :: z_colch4(klev) , z_colco2(klev)
    REAL(wp) :: z_colh2o(klev) , z_colo3(klev)
    REAL(wp) :: z_coln2o(klev) , z_colo2(klev)
    REAL(wp) :: z_forfac(klev) , z_forfrac(klev)
    REAL(wp) :: z_selffrac(klev), z_selffac(klev)
    REAL(wp) :: z_fac00(klev)  , z_fac01(klev)
    REAL(wp) :: z_fac11(klev)  , z_fac10(klev)
    REAL(wp) :: zfrcl(klev)
    REAL(wp) :: z_oneminus   , bnd_wght(ksw)
    REAL(wp) :: ztauc(klev,ksw), ztaua(klev,ksw)
    REAL(wp) :: zasyc(klev,ksw), zasya(klev,ksw)
    REAL(wp) :: zomgc(klev,ksw), zomga(klev,ksw)

    REAL(wp) :: zbbcd(klev+1,ksw), zbbcu(klev+1,ksw)
    REAL(wp) :: zbbfd(klev+1,ksw), zbbfu(klev+1,ksw)
    REAL(wp) :: zsudu(ksw), zsuduc(ksw)

    REAL(wp) :: zpm_fl_vr(klev)
    REAL(wp) :: ztk_fl_vr(klev)
    REAL(wp) :: zcol_dry_vr(klev)
    REAL(wp) :: zwkl_vr(jpinpx,klev)
    REAL(wp) :: zflxd_sw(klev+1)
    REAL(wp) :: zflxd_sw_clr(klev+1)
    REAL(wp) :: zflxd_sw_cld(klev+1)
    REAL(wp) :: zflxu_sw(klev+1)
    REAL(wp) :: zflxu_sw_clr(klev+1)
    REAL(wp) :: zflxu_sw_cld(klev+1)

    INTEGER  :: i_laytrop, i_layswtch, i_laylow
    INTEGER  :: indfor(klev), indself(klev)
    INTEGER  :: jp(klev), jt(klev), jt1(klev)

    REAL(wp) :: zclear, zcloud, zeps, zfrcl_above
    REAL(wp) :: zalbd(ksw) , zalbp(ksw)

    REAL(wp) :: frc_vis(ksw), frc_nir(ksw), frc_par
    REAL(wp) :: zfvis, zfnir, zfpar, total
    REAL(wp) :: zflxn_vis, zflxn, zflxd_vis, zflxd_nir, zflxd_par, &
                zflxd_diff, zflxd_vis_diff, zflxd_nir_diff, zflxd_par_diff, &
                zrat_swdn

    INTEGER  :: icldatm, inflag, iceflag, i_liqflag, i_nstr
    !INTEGER(i4) :: idx, >> sylvia_20200320, Warning: unused variable
    INTEGER(i4) :: ic, jk, jsw, jb, jk1, jkp1  
    ! >> sylvia_20200311, Removing jl, icount as indices here
    REAL(wp) :: zrmu0
    REAL(wp) :: ccmax, ccran, alpha, deltaz, ccrat
    LOGICAL  :: lcomp_fractions
    ! << sylvia_20200311
    
    ! >> sylvia_20200311
    ! Bringing in local variables from srtm_setcoef.f90
    INTEGER :: i_nlayers
    INTEGER :: jp1

    REAL(wp) :: z_stpfac
    REAL(wp) :: z_fp, z_ft, z_ft1, z_water, z_scalefac
    REAL(wp) :: z_factor, z_co2reg, z_compfp

    REAL(wp) :: rkindfor, rkindself
    REAL(wp) :: z_plog,z_exptavel,z_ptaveli
    ! << sylvia_20200311

    ! calculate information needed ny the radiative transfer routine
    zeps       = 1.e-06_wp
    z_oneminus = 1.0_wp - zeps

    !++hs
    ! --- weight radiation within a band for the solar cycle ---
    ! ssi contains the solar irradiation at 1 AU distance from the sun
    ! in each band. The sum over all bands of ssi(:) is the TSI.
    bnd_wght(:) = ssi(:) / ssi_default(:)
    !--hs

    icldatm   = 1
    inflag    = 2
    iceflag   = 3
    i_liqflag = 1
    i_nstr    = 2

    !-------------------------------
    ! scatter-gather in idx
    ic = 0
    ! >> sylvia_20200311
    ! Remove the loop over nproma.
    ! The threshold value of 0.0 ensures that radiation is not calculated for night points
    IF (prmu0 > 0.0_wp) THEN
      ic = ic + 1
      !idx(ic) = jl
      zrmu0 = prmu0
    ENDIF
    ! << sylvia_20200311
    ! icount=ic icount is 1. ic is 1.

    DO jk=1,klev+1
        flxu_sw(jk)     = 0.0_wp
        flxd_sw(jk)     = 0.0_wp
        flxu_sw_clr(jk) = 0.0_wp
        flxd_sw_clr(jk) = 0.0_wp
    END DO

    IF (PRESENT(flxd_dff_sfc)) THEN
        flxd_dff_sfc = 0.0_wp
    ENDIF

    IF (PRESENT(flxd_par_sfc)) THEN
        flxd_par_sfc = 0.0_wp
    ENDIF

    IF (PRESENT(vis_frc_sfc) .AND. PRESENT(nir_dff_frc_sfc) .AND. &
        PRESENT(vis_dff_frc_sfc) .AND. PRESENT(par_dff_frc_sfc)) THEN
      lcomp_fractions = .TRUE.
      vis_frc_sfc     = 0.0_wp
      nir_dff_frc_sfc = 0.0_wp
      vis_dff_frc_sfc = 0.0_wp
      par_dff_frc_sfc = 0.0_wp
    ELSE
      lcomp_fractions = .FALSE.
    END IF

    ! >> sylvia_20200320
    ! No icount variable any longer.
    !IF (icount == 0) RETURN
    ! << sylvia_20200320

    !-------------------------------

    ! >> sylvia_20200311
    ! Remove the loop over ic here.
    DO jk = 1, klev
        !jl = idx(ic)
        zfrcl(jk)       = cld_frc_vr(jk)
        zpm_fl_vr(jk)   = pm_fl_vr(jk)
        ztk_fl_vr(jk)   = tk_fl_vr(jk)
        zcol_dry_vr(jk) = col_dry_vr(jk)
!CDIR EXPAND=jpinpx
        zwkl_vr(1:jpinpx,jk) = wkl_vr(1:jpinpx,jk)
    END DO


    ! >> sylvia_20200311
    ! Removing icount, kbdim as inputs below. Implanting the subroutine itself.
    ! Start setcoef subroutine.
    z_stpfac = 296._wp/1013._wp
    i_nlayers = klev

    ! >> sylvia_20200311
    ! Removing the loop over ic
    i_layswtch = 0
    i_laytrop  = 0
    i_laylow   = 0
   

    DO jk = 1, i_nlayers
      z_plog = LOG(zpm_fl_vr(jk))
      z_ptaveli = 1.0_wp/ztk_fl_vr(jk)
      z_exptavel = z_ptaveli*EXP(-1919.4_wp*z_ptaveli)/8.7604e-4_wp

      ! Find the two reference pressures on either side of the
      ! layer pressure.  Store them in JP and JP1.  Store in FP the
      ! fraction of the difference (in ln(pressure)) between these
      ! two values that the layer pressure lies.
      jp(jk) = INT(36._wp - 5._wp*(z_plog+0.04_wp))
      IF (jp(jk) < 1) THEN
        jp(jk) = 1
      ELSEIF (jp(jk) > 58) THEN
        jp(jk) = 58
      ENDIF
      jp1 = jp(jk) + 1
      z_fp = 5._wp * (preflog(jp(jk)) - z_plog)

      !   Determine, for each reference pressure (JP and JP1), which
      !   reference temperature (these are different for each
      !   reference pressure) is nearest the layer temperature but does
      !   not exceed it.  Store these indices in JT and JT1, resp.
      !   Store in FT (resp. FT1) the fraction of the way between JT
      !   (JT1) and the next highest reference temperature that the
      !   layer temperature falls.
      jt(jk) = INT(3._wp + (ztk_fl_vr(jk)-tref(jp(jk)))/15._wp)
      IF (jt(jk) < 1) THEN
        jt(jk) = 1
      ELSEIF (jt(jk) > 4) THEN
        jt(jk) = 4
      ENDIF
      z_ft = ((ztk_fl_vr(jk)-tref(jp(jk)))/15._wp) - REAL(jt(jk)-3,wp)
      jt1(jk) = INT(3._wp + (ztk_fl_vr(jk)-tref(jp1))/15._wp)
      IF (jt1(jk) < 1) THEN
        jt1(jk) = 1
      ELSEIF (jt1(jk) > 4) THEN
        jt1(jk) = 4
      ENDIF
      z_ft1 = ((ztk_fl_vr(jk)-tref(jp1))/15._wp) - REAL(jt1(jk)-3,wp)

      z_water = zwkl_vr(1,jk)/zcol_dry_vr(jk)
      z_scalefac = zpm_fl_vr(jk) * z_stpfac * z_ptaveli

      !        If the pressure is less than ~100mb, perform a different
      !        set of species interpolations.

      IF (z_plog <= 4.56_wp) go to 5300

      i_laytrop =  i_laytrop + 1
      IF (z_plog >= 6.62_wp) i_laylow = i_laylow + 1

      !        Set up factors needed to separately include the water vapor
      !        foreign-continuum in the calculation of absorption coefficient.

      z_forfac(jk) = z_scalefac / (1._wp+z_water)
      z_factor = (332.0_wp-ztk_fl_vr(jk))/36.0_wp
      rkindfor = MIN(2.0_wp, MAX(1.0_wp, AINT(z_factor)))
      indfor(jk) = INT(rkindfor)
      z_forfrac(jk) = z_factor - rkindfor

      !        Set up factors needed to separately include the water vapor
      !        self-continuum in the calculation of absorption coefficient.

      z_selffac(jk) = z_water * z_forfac(jk)
      z_factor = (ztk_fl_vr(jk)-188.0_wp)/7.2_wp

      rkindself = MIN(9.0_wp, MAX(1.0_wp, AINT(z_factor)-7.0_wp))
      indself(jk) = INT(rkindself)
      z_selffrac(jk) = z_factor - rkindself - 7.0_wp

      !        Calculate needed column amounts.
      z_colh2o(jk) = 1.e-20_wp * zwkl_vr(1,jk)
      z_colco2(jk) = 1.e-20_wp * zwkl_vr(2,jk)
      z_colo3(jk)  = 1.e-20_wp * zwkl_vr(3,jk)
      z_coln2o(jk) = 1.e-20_wp * zwkl_vr(4,jk)
      z_colch4(jk) = 1.e-20_wp * zwkl_vr(6,jk)
      z_colo2(jk)  = 1.e-20_wp * zwkl_vr(7,jk)
      z_colmol(jk) = 1.e-20_wp * zcol_dry_vr(jk) + z_colh2o(jk)
      IF (z_colco2(jk) == 0._wp) z_colco2(jk) = 1.e-32_wp * zcol_dry_vr(jk)
      IF (z_coln2o(jk) == 0._wp) z_coln2o(jk) = 1.e-32_wp * zcol_dry_vr(jk)
      IF (z_colch4(jk) == 0._wp) z_colch4(jk) = 1.e-32_wp * zcol_dry_vr(jk)
      IF (z_colo2(jk)  == 0._wp) z_colo2(jk)  = 1.e-32_wp * zcol_dry_vr(jk)
      !        Using E = 1334.2 cm-1.
      z_co2reg = 3.55e-24_wp * zcol_dry_vr(jk)
      z_co2mult(jk)= (z_colco2(jk) - z_co2reg) * &
        & 272.63_wp * z_exptavel

      go to 5400

      !        Above LAYTROP.
5300  CONTINUE

        !        Set up factors needed to separately include the water vapor
        !        foreign-continuum in the calculation of absorption coefficient.

        z_forfac(jk) = z_scalefac / (1._wp+z_water)
        z_factor = (ztk_fl_vr(jk)-188.0_wp)/36.0_wp
        indfor(jk) = 3
        z_forfrac(jk) = z_factor - 1.0_wp

        !        Calculate needed column amounts.

        z_colh2o(jk) = 1.e-20_wp * zwkl_vr(1,jk)
        z_colco2(jk) = 1.e-20_wp * zwkl_vr(2,jk)
        z_colo3(jk)  = 1.e-20_wp * zwkl_vr(3,jk)
        z_coln2o(jk) = 1.e-20_wp * zwkl_vr(4,jk)
        z_colch4(jk) = 1.e-20_wp * zwkl_vr(6,jk)
        z_colo2(jk)  = 1.e-20_wp * zwkl_vr(7,jk)
        z_colmol(jk) = 1.e-20_wp * zcol_dry_vr(jk) + z_colh2o(jk)
        IF (z_colco2(jk) == 0._wp) z_colco2(jk) = 1.e-32_wp * zcol_dry_vr(jk)
        IF (z_coln2o(jk) == 0._wp) z_coln2o(jk) = 1.e-32_wp * zcol_dry_vr(jk)
        IF (z_colch4(jk) == 0._wp) z_colch4(jk) = 1.e-32_wp * zcol_dry_vr(jk)
        IF (z_colo2(jk)  == 0._wp) z_colo2(jk)  = 1.e-32_wp * zcol_dry_vr(jk)

        z_co2reg = 3.55e-24_wp * zcol_dry_vr(jk)
        z_co2mult(jk)= (z_colco2(jk) - z_co2reg) * &
          & 272.63_wp * z_exptavel

        z_selffac(jk) =0.0_wp
        z_selffrac(jk)=0.0_wp
        indself(jk) = 0

5400    CONTINUE

        !        We have now isolated the layer ln pressure and temperature,
        !        between two reference pressures and two reference temperatures
        !        (for each reference pressure).  We multiply the pressure
        !        fraction FP with the appropriate temperature fractions to get
        !        the factors that will be needed for the interpolation that yields
        !        the optical depths (performed in routines TAUGBn for band n).

        z_compfp = 1._wp - z_fp
        z_fac10(jk) = z_compfp * z_ft
        z_fac00(jk) = z_compfp * (1._wp - z_ft)
        z_fac11(jk) = z_fp * z_ft1
        z_fac01(jk) = z_fp * (1._wp - z_ft1)
    ENDDO

    ! >> sylvia_20200311
    ! Changing the below from finish to print.
    ! Remove the MAVAL operator and 1:count indices for i_laytrop.
    IF ( i_laytrop == klev ) THEN
      print*,('mo_srtm:srtm_setcoef: Uppermost model layer too low for RRTM.')
    ENDIF
    ! << sylvia_20200311
    ! End of setcoef subroutine.
    ! << sylvia_20200311 

    !- call the radiation transfer routine

    ! surface albedo, direct/parallel beam (p) and diffuse (d)
    DO jsw = 1, ksw
      frc_vis(jsw) = MAX(0.0_wp, MIN(1.0_wp, &
        (wavenum2(jsw+jpb1-1) - nir_vis_boundary) / delwave(jsw+jpb1-1)))
      frc_nir(jsw) = 1.0_wp - frc_vis(jsw)
!IBM* NOVECTOR
!IBM* ASSERT(NODEPS)
      ! >> sylvia_20200311
      ! Removing the loop over ic here.
      !jl = idx(ic)
      IF (nir_vis_boundary - wavenum1(jsw+jpb1-1) < 0.0_wp) THEN
        zalbd(jsw) = alb_vis_dif
        zalbp(jsw) = alb_vis_dir
      ELSEIF (nir_vis_boundary - wavenum2(jsw+jpb1-1) > 0.0_wp) THEN
        zalbd(jsw) = alb_nir_dif
        zalbp(jsw) = alb_nir_dir
      ELSE
        zalbd(jsw) = frc_vis(jsw) * alb_vis_dif &
          &           + frc_nir(jsw) * alb_nir_dif
        zalbp(jsw) = frc_vis(jsw) * alb_vis_dir &
          &           + frc_nir(jsw) * alb_nir_dir
      ENDIF
      
      ! optical properties of clouds and aerosols
      DO jk = 1, klev
!IBM* ASSERT(NODEPS)
        ! Removing the loop over ic here.
        !jl = idx(ic)
        ztauc(jk,jsw) = cld_tau_sw_vr(jsw,jk)
        zasyc(jk,jsw) = cld_cg_sw_vr(jsw,jk)
        zomgc(jk,jsw) = cld_piz_sw_vr(jsw,jk)
        ! hs: The order of indices in the aerosol optical properties was wrong. However,
        ! it should be checked if a fix in the interface would be more appropriate.
        !           ztaua(jk,jsw) = aer_tau_sw_vr(jsw,jk)
        !           zasya(jk,jsw) = aer_cg_sw_vr(jsw,jk)
        !           zomga(jk,jsw) = aer_piz_sw_vr(jsw,jk)
        ztaua(jk,jsw) = aer_tau_sw_vr(jk,jsw)
        zasya(jk,jsw) = aer_cg_sw_vr(jk,jsw)
        zomga(jk,jsw) = aer_piz_sw_vr(jk,jsw)
      ENDDO
      ! << sylvia_20200311
    ENDDO

    DO jsw=1,ksw
      DO jk=1,klev+1
!IBM* ASSERT(NODEPS)
        zbbcu(jk,jsw)=0.0_wp
        zbbcd(jk,jsw)=0.0_wp
        zbbfu(jk,jsw)=0.0_wp
        zbbfd(jk,jsw)=0.0_wp
        ! >> sylvia_20200311
        ! Removing the loop over ic here.
      ENDDO
    ENDDO


    zsudu  = 0.0_wp
    zsuduc = 0.0_wp

    ! >> sylvia_20200311
    ! Removing icount, kbdim as an input below.
    CALL srtm_spcvrt                                                       &
      !   input
      & ( klev      , ksw         ,                                        &
      &  z_oneminus,                                                       &
      &   zalbd     , zalbp       , zfrcl      ,                           &
      &   ztauc     , zasyc       , zomgc      ,                           &
      &   ztaua     , zasya       , zomga      ,                           &
      &   zrmu0     , i_laytrop   ,                                        &
      &   z_colch4  , z_colco2    , z_colh2o   ,                           &
      &   z_colmol  , z_colo2     , z_colo3    ,                           &
      &   z_forfac  , z_forfrac   , indfor     ,                           &
      &   z_selffac , z_selffrac  , indself    ,                           &
      &   z_fac00   , z_fac01     , z_fac10    , z_fac11    ,              &
      &   jp        , jt          , jt1        ,                           &
      &   zbbfd     , zbbfu       , zbbcd      , zbbcu      ,              &
      &   zsudu     , zsuduc      )

    ! >> sylvia_20200311
    ! Removing the loop over ic here.
    DO jk=1,klev+1
        zflxu_sw_cld(jk) = 0.0_wp
        zflxd_sw_cld(jk) = 0.0_wp
        zflxu_sw_clr(jk) = 0.0_wp
        zflxd_sw_clr(jk) = 0.0_wp
    END DO

    ! >> sylvia_20200311
    ! Remove (:) for zflxd_diff and zflxd_par.
    IF (PRESENT(flxd_dff_sfc)) THEN ! compute diffuse parts of surface radiation
      zflxd_diff = 0._wp
    ENDIF

    IF (PRESENT(flxd_par_sfc)) THEN ! compute photosynthetically active parts of surface radiation
      zflxd_par  = 0._wp
    ENDIF

    IF (lcomp_fractions) THEN
      zflxn_vis      = 0._wp
      zflxn          = 0._wp
      zflxd_vis      = 0._wp
      zflxd_nir      = 0._wp
      zflxd_par      = 0._wp
      zflxd_vis_diff = 0._wp
      zflxd_nir_diff = 0._wp
      zflxd_par_diff = 0._wp
    END IF

    ! >> sylvia_20200311
    ! Removing the loop over ic here.
    DO jb = 1,ksw
      ! sum up fluxes over all bands
      DO jk = 1, klev+1
!IBM* ASSERT(NODEPS)
          zflxu_sw_clr(jk) = zflxu_sw_clr(jk) + bnd_wght(jb)*zbbcu(jk,jb)
          zflxd_sw_clr(jk) = zflxd_sw_clr(jk) + bnd_wght(jb)*zbbcd(jk,jb)
          zflxu_sw_cld(jk) = zflxu_sw_cld(jk) + bnd_wght(jb)*zbbfu(jk,jb)
          zflxd_sw_cld(jk) = zflxd_sw_cld(jk) + bnd_wght(jb)*zbbfd(jk,jb)
      ENDDO
    ENDDO


    !
    ! --- overlap computation
    !
!IBM* ASSERT(NODEPS)
    ! >> sylvia_20200311
    ! Removing the loop over ic here.
    zclear      = 1.0_wp
    zcloud      = 0.0_wp
    zfrcl_above = 0.0_wp

    ! >> sylvia_20200311
    ! icld_overlap set to 2 by default according to namelists/mo_radiation_nml
    ! generalized maximum-random overlap (Hogan, Illingworth, 2000)
    ! Remove loop over ic.
    zrat_swdn = zflxd_sw_cld(klev+1)/zflxd_sw_clr(klev+1)

    DO jk = klev, 1, -1
      jkp1 = MIN(jk+1,klev)
      jk1 = klev+2-jk

      ! Remove loop over ic.
      ! reduction factor for thin cirrus clouds lying above optically thicker water clouds
      ccrat = MIN(1._wp,8._wp*MAX(5.e-5_wp,1._wp-zflxd_sw_cld(jk1)/zflxd_sw_clr(jk1))/&
                               MAX(1.e-3_wp,1._wp-zrat_swdn) )
      ccrat = MAX(ccrat,MIN(1._wp,0.1_wp*(ztk_fl_vr(jkp1)-233.15_wp))) ! no reduction for T > -30C
      ccmax = MAX( ccrat*zfrcl(jk),  zcloud )
      ccran =      ccrat*zfrcl(jk) + zcloud - zfrcl(jk) * zcloud

      ! layer thickness [m] between level jk and next upper (!) level jk+1
      deltaz = (zpm_fl_vr(jk)-zpm_fl_vr(jkp1))/(zpm_fl_vr(jkp1)+zpm_fl_vr(jk)) * &
               (ztk_fl_vr(jkp1)+ztk_fl_vr(jk))*rd*rgrav

      alpha  = MIN(EXP(-deltaz/zdecorr), zfrcl(jkp1)/MAX(zeps,zfrcl(jk)) )

      zcloud = alpha * ccmax + (1-alpha) * ccran
      zclear = 1.0_wp - zcloud
    ENDDO

    ! >> sylvia_20200311
    ! Remove loop over ic.
    DO jk = 1, klev+1
        zflxu_sw(jk) = zcloud*zflxu_sw_cld(jk) + zclear*zflxu_sw_clr(jk)
        zflxd_sw(jk) = zcloud*zflxd_sw_cld(jk) + zclear*zflxd_sw_clr(jk)
    ENDDO

    DO jb = 1,ksw

      IF (jb == 9) THEN
        frc_par = 0.533725_wp
      ELSE IF (jb == 10) THEN
        frc_par = 1.0_wp
      ELSE IF (jb == 11) THEN
        frc_par = 0.550164_wp
      ELSE
        frc_par = 0._wp
      ENDIF

      ! >> sylvia_20200311
      ! Remove loop over ic.      
      IF (PRESENT(flxd_dff_sfc)) THEN ! compute diffuse parts of surface radiation
          zflxd_diff = zflxd_diff + bnd_wght(jb)*( zcloud*(zbbfd(klev+1,jb)-zsudu(jb))  &
            &                                            + zclear*(zbbcd(klev+1,jb)-zsuduc(jb)) )
      ENDIF

      IF (PRESENT(flxd_par_sfc)) THEN ! compute photosynthetically active parts of surface radiation
          zflxd_par  = zflxd_par  + frc_par*bnd_wght(jb)*( zcloud*zbbfd(klev+1,jb) &
            &                                                    + zclear*zbbcd(klev+1,jb) )
      ENDIF

      IF (lcomp_fractions) THEN
        ! VIS, NIR and PAR fractions of bands
        zfvis = bnd_wght(jb)*frc_vis(jb)
        zfnir = bnd_wght(jb)*frc_nir(jb)
        zfpar = bnd_wght(jb)*frc_par

        zflxn_vis    = zflxn_vis + zfvis*( &
          &                 zcloud*zbbfd(klev+1,jb)     &
          &               + zclear*zbbcd(klev+1,jb)     &
          &               - zcloud*zbbfu(klev+1,jb)     &
          &               - zclear*zbbcu(klev+1,jb)     )
        zflxn = zflxd_sw(klev+1) - zflxu_sw(klev+1)
        zflxd_vis  = zflxd_vis + zfvis*(               &
          &                 zcloud*zbbfd(klev+1,jb)     &
          &               + zclear*zbbcd(klev+1,jb)     )
        zflxd_nir  = zflxd_nir + zfnir*(               &
          &                 zcloud*zbbfd(klev+1,jb)     &
          &               + zclear*zbbcd(klev+1,jb)     )
        zflxd_par  = zflxd_par + zfpar*(               &
          &                 zcloud*zbbfd(klev+1,jb)     &
          &               + zclear*zbbcd(klev+1,jb)     )
        zflxd_vis_diff = zflxd_vis_diff + zfvis*(      &
          &                 zcloud*zbbfd(klev+1,jb)     &
          &               + zclear*zbbcd(klev+1,jb)     &
          &               - zcloud*zsudu(jb)            &
          &               - zclear*zsuduc(jb))
        zflxd_nir_diff = zflxd_nir_diff + zfnir*(      &
          &                 zcloud*zbbfd(klev+1,jb)     &
          &               + zclear*zbbcd(klev+1,jb)     &
          &               - zcloud*zsudu(jb)            &
          &               - zclear*zsuduc(jb))
        zflxd_par_diff = zflxd_par_diff + zfpar*(      &
          &                 zcloud*zbbfd(klev+1,jb)     &
          &               + zclear*zbbcd(klev+1,jb)     &
          &               - zcloud*zsudu(jb)            &
          &               - zclear*zsuduc(jb))
      END IF
    END DO ! jb

      ! >> sylvia_20200311
      ! Remove loop over ic. 
      DO jk=1,klev+1
!IBM* ASSERT(NODEPS)
!CDIR NODEP,VOVERTAKE,VOB
          flxu_sw(jk)     = zflxu_sw(jk)
          flxd_sw(jk)     = zflxd_sw(jk)
          flxu_sw_clr(jk) = zflxu_sw_clr(jk)
          flxd_sw_clr(jk) = zflxd_sw_clr(jk)
      ENDDO

      ! >> sylvia_20200311
      ! Remove loop over ic. Remove jl = idx and jl indices on both fluxes below.
      IF (PRESENT(flxd_dff_sfc)) THEN ! compute diffuse parts of surface radiation
          flxd_dff_sfc = zflxd_diff
      ENDIF

      IF (PRESENT(flxd_par_sfc)) THEN ! compute photosynthetically active parts of surface radiation
          flxd_par_sfc = zflxd_par
      ENDIF

      ! >> sylvia_20200311
      ! Removing jl and ic indices here.
      IF (lcomp_fractions) THEN
          total = zflxn + zeps
          vis_frc_sfc = zflxn_vis / total
          total = zflxd_nir + zeps
          nir_dff_frc_sfc = zflxd_nir_diff / total
          total = zflxd_vis + zeps
          vis_dff_frc_sfc = zflxd_vis_diff / total
          total = zflxd_par + zeps
          par_dff_frc_sfc = zflxd_par_diff / total
      END IF


  END SUBROUTINE srtm_srtm_224gp
  
  
  SUBROUTINE srtm_spcvrt                                               &
    !   input
    & ( klev      , ksw,                                     &
    &   poneminus,                                                     &
    &   palbd   , palbp    , zfrcl     ,                               &
    &   ptauc   , pasyc    , pomgc     ,                               &
    &   ptaua   , pasya    , pomga     ,                               &
    &   zrmu0   , klaytrop ,                                           &
    &   pcolch4 , pcolco2  , pcolh2o   , pcolmol  ,                    &
    &   pcolo2  , pcolo3   ,                                           &
    &   pforfac , pforfrac , kindfor   ,                               &
    &   pselffac, pselffrac, kindself  ,                               &
    &   pfac00  , pfac01   , pfac10    , pfac11   ,                    &
    &   kjp     , kjt      , kjt1      ,                               &
    !   output
    &   pbbfd   , pbbfu    , pbbcd     , pbbcu    ,                    &
    &   psudu   , psuduc )

    !**** *SRTM_SPCVRT* - SPECTRAL LOOP TO COMPUTE THE SHORTWAVE RADIATION FLUXES.

    !     PURPOSE.
    !     --------

    !          THIS ROUTINE COMPUTES THE TWO-STREAM METHOD OF BARKER

    !**   INTERFACE.
    !     ----------

    !          *SRTM_SPCVRT* IS CALLED FROM *SRTM_SRTM_224GP*

    !        IMPLICIT ARGUMENTS :
    !        --------------------

    !     ==== INPUTS ===
    !     ==== OUTPUTS ===

    !     METHOD.
    !     -------

    !     EXTERNALS.
    !     ----------

    !          *SWVRTQDR*

    !     REFERENCE.
    !     ----------

    !        SEE RADIATION'S PART OF THE ECMWF RESEARCH DEPARTMENT
    !        DOCUMENTATION
    !     AUTHOR.
    !     -------
    !        from Howard Barker
    !        JEAN-JACQUES MORCRETTE  *ECMWF*

    !     MODIFICATIONS.
    !     --------------
    !        ORIGINAL : 03-02-27
    !        M.Hamrud      01-Oct-2003 CY28 Cleaning
    !        JJMorcrette 20070614 bug-fix for solar duration
    !        D.Salmond  31-Oct-2007 Vector version in the style of RRTM from Meteo France & NEC
    !        M.Puetz    20-Apr-2010 manual gather/scatter for prmu0 > 0, index passed to subroutines
    !     ------------------------------------------------------------------

    !*       0.1   ARGUMENTS
    !              ---------

    !-- Input arguments
    
    ! >> sylvia_20200310
    ! Set up the wp type that is often used below.
    ! Pull the nbndlw and ngptlw constants from mo_lrtm_par.f90
    INTEGER, PARAMETER :: pd =  12
    INTEGER, PARAMETER :: rdd = 307
    INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(pd,rdd) !< double precision
    INTEGER, PARAMETER :: wp = dp
    INTEGER, PARAMETER :: pi4 =  9
    INTEGER, PARAMETER :: i4 = SELECTED_INT_KIND(pi4)   !< at least 4 byte integer
    ! << sylvia_20200311
    
    ! << sylvia_20200311
    ! Removing icount here.
    INTEGER,  INTENT(in)    :: ksw
    INTEGER,  INTENT(in)    :: klev
    REAL(wp), INTENT(in)    :: poneminus
    REAL(wp), INTENT(in)    :: palbd(ksw)
    REAL(wp), INTENT(in)    :: palbp(ksw)
    REAL(wp), INTENT(in)    :: zfrcl(klev)     ! bottom to top
    REAL(wp), INTENT(in)    :: ptauc(klev,ksw) ! bottom to top
    REAL(wp), INTENT(in)    :: pasyc(klev,ksw) ! bottom to top
    REAL(wp), INTENT(in)    :: pomgc(klev,ksw) ! bottom to top
    REAL(wp), INTENT(in)    :: ptaua(klev,ksw) ! bottom to top
    REAL(wp), INTENT(in)    :: pasya(klev,ksw) ! bottom to top
    REAL(wp), INTENT(in)    :: pomga(klev,ksw) ! bottom to top
    REAL(wp), INTENT(in)    :: zrmu0
    INTEGER,  INTENT(in)    :: klaytrop
    REAL(wp), INTENT(in)    :: pcolch4(klev)
    REAL(wp), INTENT(in)    :: pcolco2(klev)
    REAL(wp), INTENT(in)    :: pcolh2o(klev)
    REAL(wp), INTENT(in)    :: pcolmol(klev)
    REAL(wp), INTENT(in)    :: pcolo2(klev)
    REAL(wp), INTENT(in)    :: pcolo3(klev)
    REAL(wp), INTENT(in)    :: pforfac(klev)
    REAL(wp), INTENT(in)    :: pforfrac(klev)
    INTEGER,  INTENT(in)    :: kindfor(klev)
    REAL(wp), INTENT(in)    :: pselffac(klev)
    REAL(wp), INTENT(in)    :: pselffrac(klev)
    INTEGER,  INTENT(in)    :: kindself(klev)
    REAL(wp), INTENT(in)    :: pfac00(klev)
    REAL(wp), INTENT(in)    :: pfac01(klev)
    REAL(wp), INTENT(in)    :: pfac10(klev)
    REAL(wp), INTENT(in)    :: pfac11(klev)
    INTEGER,  INTENT(in)    :: kjp(klev)
    INTEGER,  INTENT(in)    :: kjt(klev)
    INTEGER,  INTENT(in)    :: kjt1(klev)

    !-- Input and output arguments

    REAL(wp) ,INTENT(inout) :: pbbfd(klev+1,jpb1:jpb2)
    REAL(wp) ,INTENT(inout) :: pbbfu(klev+1,jpb1:jpb2)
    REAL(wp) ,INTENT(inout) :: pbbcd(klev+1,jpb1:jpb2)
    REAL(wp) ,INTENT(inout) :: pbbcu(klev+1,jpb1:jpb2)
    REAL(wp) ,INTENT(inout) :: psudu(jpb1:jpb2)
    REAL(wp) ,INTENT(inout) :: psuduc(jpb1:jpb2)

    !-- Local variables

    LOGICAL :: llrtchk(klev)

    REAL(wp) ::  zclear  , zcloud
    REAL(wp) ::                                               &
      &   zdbt(klev+1)                                        &
      & , zgcc(klev)   , zgco(klev)                           &
      & , zomcc(klev)  , zomco(klev)                          &
      & , zrdnd(klev+1), zrdndc(klev+1)                       &
      & , zref(klev+1) , zrefc(klev+1) , zrefo(klev+1)  &
      & , zrefd(klev+1), zrefdc(klev+1), zrefdo(klev+1) &
      & , zrup(klev+1) , zrupd(klev+1)                        &
      & , zrupc(klev+1), zrupdc(klev+1)                       &
      & , ztauc(klev)  , ztauo(klev)                          &
      & , ztdbt(klev+1)                                       &
      & , ztra(klev+1) , ztrac(klev+1) , ztrao(klev+1)        &
      & , ztrad(klev+1), ztradc(klev+1), ztrado(klev+1)
    REAL(wp) ::                                               &
      & zdbtc(klev+1), ztdbtc(klev+1)  , zincflx              &
      & , zincf14(14)  , zinctot

    INTEGER :: ib1, ib2, ibm, igt, ikl, jb, jg, jk, i_kmodts
    REAL(wp) :: zdbtmc, zdbtmo, zf, zincflux, zwf

    !-- Output of SRTM_TAUMOLn routines

    REAL(wp) :: ztaug(klev,16), ztaur(klev,16), zsflxzen(16)

    !-- Output of SRTM_VRTQDR routine
    REAL(wp) ::                                            &
      &   zcd(klev+1), zcu(klev+1) &
      & , zfd(klev+1), zfu(klev+1)

    REAL(wp) :: zrmu0i
    ! >> sylvia_20200311
    ! Removing ic integer here.

    ! >> sylvia_20200320
    ! Adding inputs back in from srt_vrtqdr subroutine without INTENT(in)
    ! Taking out pdbt, pfd, pfu, prdnd, pref, prefd, prup, prupd, ptdbt,
    ! ptra, ptrad
    REAL(wp) :: ztdn(klev+1)
    INTEGER(i4) :: ikp, ikx   !, ic >> sylvia_20200320, Unused variable.
    REAL(wp) :: zreflect
    ! << sylvia_20200320

    !     ------------------------------------------------------------------

    !-- Two-stream model 1: Eddington, 2: PIFM, Zdunkowski et al., 3: discret ordinates
    ! KMODTS is set in SWREFTRA
    !NDBUG=4
    ibm = -1
    igt = -1
    ib1=jpb1
    ib2=jpb2

!IBM* ASSERT(NODEPS)
    ! >> sylvia_20200311
    ! Remove the loop over ic.
    zincflux=0.0_wp
    zinctot =0.0_wp
    zrmu0i  = 1._wp / zrmu0
    ! << sylvia_20203011

    jb = ib1 - 1
    DO jb = ib1, ib2
!IBM* ASSERT(NODEPS)
      
      ibm = jb-15
      igt = ngc(ibm)
      zincf14(ibm)=0.0_wp
      
      !-- for each band, computes the gaseous and Rayleigh optical thickness
      !  for all g-points within the band
      ! >> sylvia_20200311
      ! Removing icount, kbdim as an input for all of the below.
      IF (jb == 16) THEN
        CALL srtm_taumol16                                                          &
          & ( klev    ,                                                             &
          &   pfac00  , pfac01  , pfac10   , pfac11   ,                             &
          &   kjp     , kjt     , kjt1     , poneminus,                             &
          &   pcolh2o , pcolch4 , pcolmol  ,                                        &
          &   klaytrop, pselffac, pselffrac, kindself , pforfac, pforfrac, kindfor, &
          &   zsflxzen, ztaug   , ztaur                                      &
          & )

      ELSEIF (jb == 17) THEN
        CALL srtm_taumol17                                                          &
          & ( klev    ,                                                             &
          &   pfac00  , pfac01  , pfac10   , pfac11   ,                             &
          &   kjp     , kjt     , kjt1     , poneminus,                             &
          &   pcolh2o , pcolco2 , pcolmol  ,                                        &
          &   klaytrop, pselffac, pselffrac, kindself , pforfac, pforfrac, kindfor, &
          &   zsflxzen, ztaug   , ztaur                                     &
          & )

      ELSEIF (jb == 18) THEN
        CALL srtm_taumol18                                                          &
          & ( klev    ,                                                   &
          &   pfac00  , pfac01  , pfac10   , pfac11   ,                             &
          &   kjp     , kjt     , kjt1     , poneminus,                             &
          &   pcolh2o , pcolch4 , pcolmol  ,                                        &
          &   klaytrop, pselffac, pselffrac, kindself , pforfac, pforfrac, kindfor, &
          &   zsflxzen, ztaug   , ztaur                                     &
          & )

      ELSEIF (jb == 19) THEN
        CALL srtm_taumol19                                                          &
          & ( klev    ,                                                   &
          &   pfac00  , pfac01  , pfac10   , pfac11   ,                             &
          &   kjp     , kjt     , kjt1     , poneminus,                             &
          &   pcolh2o , pcolco2 , pcolmol  ,                                        &
          &   klaytrop, pselffac, pselffrac, kindself , pforfac, pforfrac, kindfor, &
          &   zsflxzen, ztaug   , ztaur                                  &
          & )

      ELSEIF (jb == 20) THEN
        CALL srtm_taumol20                                                          &
          & ( klev    ,                                                   &
          &   pfac00  , pfac01  , pfac10   , pfac11   ,                             &
          &   kjp     , kjt     , kjt1     ,                                        &
          &   pcolh2o , pcolch4 , pcolmol  ,                                        &
          &   klaytrop, pselffac, pselffrac, kindself , pforfac, pforfrac, kindfor, &
          &   zsflxzen, ztaug   , ztaur                                      &
          & )

      ELSEIF (jb == 21) THEN
        CALL srtm_taumol21                                                          &
          & ( klev    ,                                                   &
          &   pfac00  , pfac01  , pfac10   , pfac11   ,                             &
          &   kjp     , kjt     , kjt1     , poneminus,                             &
          &   pcolh2o , pcolco2 , pcolmol  ,                                        &
          &   klaytrop, pselffac, pselffrac, kindself , pforfac, pforfrac, kindfor, &
          &   zsflxzen, ztaug   , ztaur                                      &
          & )

      ELSEIF (jb == 22) THEN
        CALL srtm_taumol22                                                          &
          & ( klev    ,                                                   &
          &   pfac00  , pfac01  , pfac10   , pfac11   ,                             &
          &   kjp     , kjt     , kjt1     , poneminus,                             &
          &   pcolh2o , pcolmol , pcolo2   ,                                        &
          &   klaytrop, pselffac, pselffrac, kindself , pforfac, pforfrac, kindfor, &
          &   zsflxzen, ztaug   , ztaur                                     &
          & )

      ELSEIF (jb == 23) THEN
        CALL srtm_taumol23                                                          &
          & ( klev    ,                                                   &
          &   pfac00  , pfac01  , pfac10   , pfac11   ,                             &
          &   kjp     , kjt     , kjt1     ,                                        &
          &   pcolh2o , pcolmol ,                                                   &
          &   klaytrop, pselffac, pselffrac, kindself , pforfac, pforfrac, kindfor, &
          &   zsflxzen, ztaug   , ztaur                                     &
          & )

      ELSEIF (jb == 24) THEN
        CALL srtm_taumol24                                                          &
          & ( klev     ,                                        &
          &   pfac00  , pfac01  , pfac10   , pfac11   ,                             &
          &   kjp     , kjt     , kjt1     , poneminus,                             &
          &   pcolh2o , pcolmol , pcolo2   , pcolo3   ,                             &
          &   klaytrop, pselffac, pselffrac, kindself , pforfac, pforfrac, kindfor, &
          &   zsflxzen, ztaug   , ztaur                                     &
          & )

      ELSEIF (jb == 25) THEN
        !--- visible 16000-22650 cm-1   0.4415 - 0.6250 um
        CALL srtm_taumol25                              &
          & ( klev     ,            &
          &   pfac00  , pfac01  , pfac10   , pfac11   , &
          &   kjp     , kjt     , kjt1     ,            &
          &   pcolh2o , pcolmol , pcolo3   ,            &
          &   klaytrop,                                 &
          &   zsflxzen, ztaug   , ztaur          &
          & )

      ELSEIF (jb == 26) THEN
        !--- UV-A 22650-29000 cm-1   0.3448 - 0.4415 um
        CALL srtm_taumol26                              &
          & ( klev  ,       &
          &   pcolmol ,                                 &
          &   klaytrop,                                 &
          &   zsflxzen, ztaug   , ztaur          &
          & )

      ELSEIF (jb == 27) THEN
        !--- UV-B 29000-38000 cm-1   0.2632 - 0.3448 um
        CALL srtm_taumol27                              &
          & ( klev     ,    &
          &   pfac00  , pfac01  , pfac10   , pfac11   , &
          &   kjp     , kjt     , kjt1     ,            &
          &   pcolmol , pcolo3  ,                       &
          &   klaytrop,                                 &
          &   zsflxzen, ztaug   , ztaur                 &
          & )

      ELSEIF (jb == 28) THEN
        !--- UV-C 38000-50000 cm-1   0.2000 - 0.2632 um
        CALL srtm_taumol28 &
          & ( klev     ,   &
          &   pfac00  , pfac01  , pfac10   , pfac11   , &
          &   kjp     , kjt     , kjt1     , poneminus, &
          &   pcolmol , pcolo2  , pcolo3   ,            &
          &   klaytrop,                                 &
          &   zsflxzen, ztaug   , ztaur                 &
          & )

      ELSEIF (jb == 29) THEN
        CALL srtm_taumol29                                                          &
          & ( klev     ,                                        &
          &   pfac00  , pfac01  , pfac10   , pfac11   ,                             &
          &   kjp     , kjt     , kjt1     ,                                        &
          &   pcolh2o , pcolco2 , pcolmol  ,                                        &
          &   klaytrop, pselffac, pselffrac, kindself , pforfac, pforfrac, kindfor, &
          &   zsflxzen, ztaug   , ztaur                                     &
          & )

      ENDIF

      DO jg=1,igt

!IBM* ASSERT(NODEPS)

          ! >> sylvia_20200311
          ! Remove the loop over ic.
          zincflx     =zsflxzen(jg)*zrmu0
          zincflux    =zincflux+zsflxzen(jg)*zrmu0
          zinctot     =zinctot+zsflxzen(jg)
          zincf14(ibm) =zincf14(ibm)+zsflxzen(jg)
            !-- CALL to compute layer reflectances and transmittances for direct
            !  and diffuse sources, first clear then cloudy.
            !   Use direct/parallel albedo for direct radiation and diffuse albedo
            !   otherwise.

            ! ZREFC  direct albedo for clear
            ! ZREFO  direct albedo for cloud
            ! ZREFDC diffuse albedo for clear
            ! ZREFDO diffuse albedo for cloud
            ! ZTRAC  direct transmittance for clear
            ! ZTRAO  direct transmittance for cloudy
            ! ZTRADC diffuse transmittance for clear
            ! ZTRADO diffuse transmittance for cloudy

            ! ZREF   direct reflectance
            ! ZREFD  diffuse reflectance
            ! ZTRA   direct transmittance
            ! ZTRAD  diffuse transmittance

            ! ZDBTC  clear direct beam transmittance
            ! ZDBTO  cloudy direct beam transmittance
            ! ZDBT   layer mean direct beam transmittance
            ! ZTDBT  total direct beam transmittance at levels

            !-- clear-sky
            !----- TOA direct beam
            ztdbtc(1)=1._wp
            !----- surface values
            zdbtc(klev+1) =0.0_wp
            ztrac(klev+1) =0.0_wp
            ztradc(klev+1)=0.0_wp
            zrefc(klev+1) =palbp(ibm)
            zrefdc(klev+1)=palbd(ibm)
            zrupc(klev+1) =palbp(ibm)
            zrupdc(klev+1)=palbd(ibm)

            !-- total sky
            !----- TOA direct beam
            ztdbt(1)=1._wp
            !----- surface values
            zdbt(klev+1) =0.0_wp
            ztra(klev+1) =0.0_wp
            ztrad(klev+1)=0.0_wp
            zref(klev+1) =palbp(ibm)
            zrefd(klev+1)=palbd(ibm)
            zrup(klev+1) =palbp(ibm)
            zrupd(klev+1)=palbd(ibm)

!IBM* NOUNROLL
        DO jk=1,klev
            ikl=klev+1-jk
            !-- NB: a two-stream calculations from top to bottom, but RRTM_SW quantities
            !       are given bottom to top (argh!)
            !       Inputs for clouds and aerosols are bottom to top as inputs

            !-- clear-sky optical parameters
            ! llrtchk(jk)=.TRUE. here always true

            !-- clear-sky optical parameters including aerosols
            ztauc(jk) = ztaur(ikl,jg) + ztaug(ikl,jg) + ptaua(ikl,ibm)
            ztauo(jk) = ztauc(jk) + ptauc(ikl,ibm)

            zomcc(jk) = ztaur(ikl,jg) + ptaua(ikl,ibm)*pomga(ikl,ibm)
            zomco(jk) = zomcc(jk) + ptauc(ikl,ibm)*pomgc(ikl,ibm)

            zgcc(jk)  = pasya(ikl,ibm)*pomga(ikl,ibm)*ptaua(ikl,ibm)
            zgco(jk)  = zgcc(jk) + ptauc(ikl,ibm)*pomgc(ikl,ibm)*pasyc(ikl,ibm)

            zgcc(jk)  = zgcc(jk)  / zomcc(jk)
            zgco(jk)  = zgco(jk)  / zomco(jk)
            zomcc(jk) = zomcc(jk) / ztauc(jk)
            zomco(jk) = zomco(jk) / ztauo(jk)

!IBM* NOVECTOR

            !-- Delta scaling for clear-sky / aerosol optical quantities
            zf=zgcc(jk)*zgcc(jk)
            zwf=zomcc(jk)*zf
            ztauc(jk)=MAX(0._wp,(1._wp-zwf)*ztauc(jk))
            zomcc(jk)=(zomcc(jk)-zwf)/(1.0_wp-zwf)
            zgcc(jk)=(zgcc(jk)-zf)/(1.0_wp-zf)

            !-- Delta scaling for cloudy quantities
            zf=zgco(jk)*zgco(jk)
            zwf=zomco(jk)*zf
            ztauo(jk)=MAX(0._wp,(1._wp-zwf)*ztauo(jk))
            zomco(jk)=(zomco(jk)-zwf)/(1._wp-zwf)
            zgco(jk)=(zgco(jk)-zf)/(1._wp-zf)
        ENDDO
        ! mpuetz: llrtchk = TRUE always in this call to reftra()
        ! >> sylvia_20200320
        ! Remove icount and kbdim as inputs to this subroutine.
        CALL srtm_reftra ( klev, i_kmodts ,&
          &   zgcc  , zrmu0, zrmu0i, ztauc , zomcc ,&
          &   zrefc  , zrefdc, ztrac, ztradc )
        ! << sylvia_20200320

        ! >> sylvia_20200310
        ! Loop over ic removed.
        DO jk = 1,klev
          ikli = klev + 1 - jk
          llrtchk(jk) = (zfrcl(ikl) > repclc)
        END DO

        CALL srtm_reftra ( klev, i_kmodts ,&
          &   zgco  , zrmu0, zrmu0i, ztauo , zomco ,&
          &   zrefo , zrefdo, ztrao, ztrado , llrtchk )

        DO jk=1,klev
          ikl=klev+1-jk
!IBM* ASSERT(NODEPS)

            !-- combine clear and cloudy contributions for total sky
            zclear   = 1.0_wp - zfrcl(ikl)
            zcloud   = zfrcl(ikl)

            zref(jk) = zclear*zrefc(jk) + zcloud*zrefo(jk)
            zrefd(jk)= zclear*zrefdc(jk)+ zcloud*zrefdo(jk)
            ztra(jk) = zclear*ztrac(jk) + zcloud*ztrao(jk)
            ztrad(jk)= zclear*ztradc(jk)+ zcloud*ztrado(jk)

            !-- direct beam transmittance

            zdbtmc     = EXP(-ztauc(jk)*zrmu0i)
            zdbtmo     = EXP(-ztauo(jk)*zrmu0i)
            zdbt(jk)   = zclear*zdbtmc+zcloud*zdbtmo
            ztdbt(jk+1)= zdbt(jk)*ztdbt(jk)

            !-- clear-sky
            zdbtc(jk)   =zdbtmc
            ztdbtc(jk+1)=zdbtc(jk)*ztdbtc(jk)
        ENDDO


        !-- vertical quadrature producing clear-sky fluxes
        ! >> sylvia_20200311
        ! Start srtm_vrtqdr in for the first time.
!IBM* ASSERT(NODEPS)
!IBM* NOVECTOR
        ! >> sylvia_20200311
        ! Remove the loop over icount
        zreflect = 1.0_wp / (1.0_wp -zrefdc(klev+1)*zrefdc(klev))
        zrupc(klev) = zrefc(klev)+(ztradc(klev)* &
             & ((ztrac(klev)-zdbtc(klev))*zrefdc(klev+1)+ &
             & zdbtc(klev)*zrefc(klev+1)))*zreflect
        zrupdc(klev)=zrefdc(klev)+ztradc(klev)* &
             & ztradc(klev)*zrefdc(klev+1)*zreflect
        ! << sylvia_20200311

        !-- pass from bottom to top

        DO jk=1,klev-1
            ikp=klev+1-jk
            ikx=ikp-1
!IBM* ASSERT(NODEPS)
!IBM* NOVECTOR
            zreflect=1.0_wp / (1.0_wp -zrupdc(ikp)*zrefdc(ikx))
            zrupc(ikx)=zrefc(ikx)+(ztradc(ikx)* &
                & ((ztrac(ikx)-zdbtc(ikx))*zrupdc(ikp)+ &
                & zdbtc(ikx)*zrupc(ikp)))*zreflect
            zrupdc(ikx)=zrefdc(ikx)+ztradc(ikx)* &
                & ztradc(ikx)*zrupdc(ikp)*zreflect
        ENDDO

        !-- upper boundary conditions

!IBM* ASSERT(NODEPS)
        ztdn(1)=1.0_wp
        zrdndc(1)=0.0_wp
        ztdn(2)=ztrac(1)
        zrdndc(2)=zrefdc(1)

        !-- pass from top to bottom
        DO jk=2,klev
            ikp=jk+1
!IBM* ASSERT(NODEPS)
!IBM* NOVECTOR
            zreflect=1.0_wp / (1.0_wp -zrefdc(jk)*zrdndc(jk))
            ztdn(ikp)=ztdbtc(jk)*ztrac(jk)+ &
                & (ztradc(jk)*((ztdn(jk)-ztdbtc(jk))+ &
                & ztdbtc(jk)*zrefc(jk)*zrdndc(jk))) * zreflect
            zrdndc(ikp)=zrefdc(jk)+ztradc(jk)*ztradc(jk) &
                & *zrdndc(jk)*zreflect
        ENDDO

        !-- up and down-welling fluxes at levels
        DO jk=1,klev+1
!IBM* ASSERT(NODEPS)
!IBM* NOVECTOR
            zreflect=1.0_wp / (1.0_wp - zrdndc(jk)*zrupdc(jk))
            zcu(jk)=(ztdbtc(jk)*zrupc(jk) + &
                & (ztdn(jk)-ztdbtc(jk))*zrupdc(jk))*zreflect
            zcd(jk)=ztdbtc(jk) + (ztdn(jk)-ztdbtc(jk)+ &
                & ztdbtc(jk)*zrupc(jk)*zrdndc(jk))*zreflect
        ENDDO
        ! End of srtm_vrtqdr for the first time.
        ! << sylvia_20200311

        !-- vertical quadrature producing cloudy fluxes
        ! Start srtm_vrtqdr for the second time.
!IBM* ASSERT(NODEPS)
!IBM* NOVECTOR
        ! >> sylvia_20200311
        ! Remove the loop over icount
        zreflect = 1.0_wp / (1.0_wp -zrefd(klev+1)*zrefd(klev))
        zrup(klev) = zref(klev)+(ztrad(klev)* &
            & ((ztra(klev)-zdbt(klev))*zrefd(klev+1)+ &
            & zdbt(klev)*zref(klev+1)))*zreflect
        zrupd(klev)=zrefd(klev)+ztrad(klev)* &
            & ztrad(klev)*zrefd(klev+1)*zreflect
        ! << sylvia_20200311

        !-- pass from bottom to top

        DO jk=1,klev-1
            ikp=klev+1-jk
            ikx=ikp-1
!IBM* ASSERT(NODEPS)
!IBM* NOVECTOR
            zreflect=1.0_wp / (1.0_wp -zrupd(ikp)*zrefd(ikx))
            zrup(ikx)=zref(ikx)+(ztrad(ikx)* &
                & ((ztra(ikx)-zdbt(ikx))*zrupd(ikp)+ &
                & zdbt(ikx)*zrup(ikp)))*zreflect
            zrupd(ikx)=zrefd(ikx)+ztrad(ikx)* &
                & ztrad(ikx)*zrupd(ikp)*zreflect
        ENDDO

        !-- upper boundary conditions

!IBM* ASSERT(NODEPS)
        ztdn(1)=1.0_wp
        zrdnd(1)=0.0_wp
        ztdn(2)=ztra(1)
        zrdnd(2)=zrefd(1)

        !-- pass from top to bottom
        DO jk=2,klev
            ikp=jk+1
!IBM* ASSERT(NODEPS)
!IBM* NOVECTOR
            zreflect=1.0_wp / (1.0_wp -zrefd(jk)*zrdnd(jk))
            ztdn(ikp)=ztdbt(jk)*ztra(jk)+ &
                & (ztrad(jk)*((ztdn(jk)-ztdbt(jk))+ &
                & ztdbt(jk)*zref(jk)*zrdnd(jk))) * zreflect
            zrdnd(ikp)=zrefd(jk)+ztrad(jk)*ztrad(jk) &
                & *zrdnd(jk)*zreflect
        ENDDO

    !-- up and down-welling fluxes at levels
        DO jk=1,klev+1
!IBM* ASSERT(NODEPS)
!IBM* NOVECTOR
            zreflect=1.0_wp / (1.0_wp - zrdnd(jk)*zrupd(jk))
            zfu(jk)=(ztdbt(jk)*zrup(jk) + &
                & (ztdn(jk)-ztdbt(jk))*zrupd(jk))*zreflect
            zfd(jk)=ztdbt(jk) + (ztdn(jk)-ztdbt(jk)+ &
                & ztdbt(jk)*zrup(jk)*zrdnd(jk))*zreflect
        ENDDO
        ! >> sylvia_20200311
        ! End srtm_vrtqdr for the second time.


        !-- up and down-welling fluxes at levels
        ! >> sylvia_20200320
        ! Assume that icount > 0.
        !IF (icount > 0) THEN
!IBM* NOUNROLL
        DO jk=1,klev+1
!IBM* ASSERT(NODEPS)
          ! >> sylvia_20200320
          ! Commenting out the ic loop. 
          !DO ic = 1,icount
          !-- accumulation of spectral fluxes
          pbbfu(jk,jb) = pbbfu(jk,jb) + zincflx*zfu(jk)
          pbbfd(jk,jb) = pbbfd(jk,jb) + zincflx*zfd(jk)
          pbbcu(jk,jb) = pbbcu(jk,jb) + zincflx*zcu(jk)
          pbbcd(jk,jb) = pbbcd(jk,jb) + zincflx*zcd(jk)
          !ENDDO
          ! >> sylvia_20200320
        ENDDO
        !END IF
!IBM* ASSERT(NODEPS)
        ! >> sylvia_20200320
        ! Commenting out the ic loop.
        !DO ic = 1,icount
        psudu(jb) = psudu(jb)  +zincflx*ztdbt (klev+1)
        psuduc(jb)= psuduc(jb) +zincflx*ztdbtc(klev+1)
        !ENDDO
        !
        ! -- direct flux
        !

      ENDDO
      !-- end loop on JG

    ENDDO
    !-- end loop on JB
    !------------------------------------------------------------------

  END SUBROUTINE srtm_spcvrt
  
 SUBROUTINE srtm_reftra &
       & ( klev  , kmodts, &
       &   pgg   , prmuz, prmuzi, ptau , pw, &
       &   pref  , prefd, ptra , ptrad , &
       &   ldrtchk &
       & )

    !**** *SRTM_REFTRA* - REFLECTIVITY AND TRANSMISSIVITY

    !     PURPOSE.
    !     --------
    !           COMPUTES THE REFLECTIVITY AND TRANSMISSIVITY OF A CLEAR OR
    !     CLOUDY LAYER USING A CHOICE OF VARIOUS APPROXIMATIONS.

    !**   INTERFACE.
    !     ----------
    !          *SRTM_REFTRA* IS CALLED BY *SRTM_SPCVRT*

    !        EXPLICIT ARGUMENTS :
    !        --------------------
    ! INPUTS
    ! ------
    !      KMODTS  = 1 EDDINGTON (JOSEPH ET AL., 1976)
    !              = 2 PIFM (ZDUNKOWSKI ET AL., 1980)
    !              = 3 DISCRETE ORDINATES (LIOU, 1973)
    !      LDRTCHK = .T. IF CLOUDY
    !              = .F. IF CLEAR-SKY
    !      PGG     = ASSYMETRY FACTOR
    !      PRMUZ   = COSINE SOLAR ZENITH ANGLE
    !      PRMUZI  = INVERSE COSINE SOLAR ZENITH ANGLE
    !      PTAU    = OPTICAL THICKNESS
    !      PW      = SINGLE SCATTERING ALBEDO

    ! OUTPUTS
    ! -------
    !      PREF    : COLLIMATED BEAM REFLECTIVITY
    !      PREFD   : DIFFUSE BEAM REFLECTIVITY
    !      PTRA    : COLLIMATED BEAM TRANSMISSIVITY
    !      PTRAD   : DIFFUSE BEAM TRANSMISSIVITY

    !     METHOD.
    !     -------
    !          STANDARD DELTA-EDDINGTON, P.I.F.M., OR D.O.M. LAYER CALCULATIONS.

    !     EXTERNALS.
    !     ----------
    !          NONE

    !     REFERENCE.
    !     ----------

    !     AUTHOR.
    !     -------
    !        JEAN-JACQUES MORCRETTE  *ECMWF*

    !     MODIFICATIONS.
    !     --------------
    !        ORIGINAL : 03-02-27
    !        M.Hamrud   01-Oct-2003      CY28 Cleaning
    !        Mike Iacono, AER, Mar 2004: bug fix
    !        D.Salmond  31-Oct-2007 Vector version in the style of RRTM from Meteo France & NEC
    !        M. Puetz   20-Apr-2010 Gather/Scatter for better pipelining on scalar cpus
    !     ------------------------------------------------------------------

    !*       0.1   ARGUMENTS
    !              ---------
    INTEGER,PARAMETER         :: icount = 1
    INTEGER,INTENT(in)        :: klev
    INTEGER,INTENT(out)       :: kmodts
    REAL(wp)   ,INTENT(in)    :: pgg(klev)
    REAL(wp)   ,INTENT(in)    :: prmuz
    REAL(wp)   ,INTENT(in)    :: prmuzi
    REAL(wp)   ,INTENT(in)    :: ptau(klev)
    REAL(wp)   ,INTENT(in)    :: pw(klev)
    REAL(wp)   ,INTENT(out)   :: pref(klev)
    REAL(wp)   ,INTENT(inout) :: prefd(klev)
    REAL(wp)   ,INTENT(inout) :: ptra(klev)
    REAL(wp)   ,INTENT(inout) :: ptrad(klev)
    LOGICAL,INTENT(in),optional :: ldrtchk(klev)
    !     ------------------------------------------------------------------
    INTEGER(i4) :: idxt, idxf, idxc, idxn
    INTEGER(i4) :: jk, ic, jc, ict, icf, icc, icn

    REAL(wp) :: zgamma1,zgamma2,zgamma3,zcrit
    REAL(wp) :: zrk, zem1, zem2, zep1, zep2
    REAL(wp) :: za, za1, za2, zemm
    REAL(wp) :: zbeta, zdend, zdenr, zdent
    REAL(wp) :: zg, zg3, zgamma4, zgt
    REAL(wp) :: zr1, zr2, zr3, zr4, zr5, zrk2, zrkg, zrm1, zrp, zrp1, zrpp
    REAL(wp) :: zsr3, zt1, zt2, zt3 ! , zt4, zt5
    REAL(wp) :: zw, zwcrit
    REAL(wp) :: ztemp, zzz, zeps, zaux
    !------------------------------------------------------------------

    zeps = 1.e-20_wp

    zsr3=SQRT(3._wp)
    zwcrit=0.9995_wp
    kmodts=2
    ict = -1

    IF (.NOT.PRESENT(ldrtchk)) THEN
      ict = icount
      icf = 0
      DO ic=1,icount
        idxt = ic
      END DO
    END IF

!PREVENT_INCONSISTENT_IFORT_FMA
    DO jk=1,klev
      IF (PRESENT(ldrtchk)) THEN
        ict = 0
        icf = 0
        DO ic=1,icount
          IF (ldrtchk(jk)) THEN
            ict = ict + 1
            idxt = ic
          ELSE
            icf = icf + 1
            idxf = ic
          END IF
        END DO
      END IF
       !-- GENERAL TWO-STREAM EXPRESSIONS

      IF (kmodts == 1) THEN
!IBM* ASSERT(NODEPS)
!CDIR NODEP,VOVERTAKE,VOB
        DO jc = 1,ict
          ic = idxt

          zw  =pw(jk)
          zg  =pgg(jk)

          zg3= 3._wp * zg
          zgamma1= (7._wp - zw * (4._wp + zg3)) * 0.25_wp
          zgamma2=-(1._wp - zw * (4._wp - zg3)) * 0.25_wp
          zgamma3= (2._wp - zg3 * prmuz ) * 0.25_wp
          zzz=(1._wp - zg)**2
          zcrit = zw*zzz - zwcrit*(zzz - (1._wp - zw)*(zg **2))
        END DO
      ELSEIF (kmodts == 2 .AND. ict == icount) THEN ! use direct addressing
        DO ic = 1,ict

          zw  =pw(jk)
          zg  =pgg(jk)

          zg3= 3._wp * zg
          zgamma1= (8._wp - zw * (5._wp + zg3)) * 0.25_wp
          zgamma2=  3._wp *(zw * (1._wp - zg )) * 0.25_wp
          zgamma3= (2._wp - zg3 * prmuz ) * 0.25_wp
          zzz=(1._wp - zg)**2
          zcrit = zw*zzz - zwcrit*(zzz - (1._wp - zw)*(zg **2))
        END DO
      ELSEIF (kmodts == 2) THEN
!IBM* ASSERT(NODEPS)
!CDIR NODEP,VOVERTAKE,VOB
!DIR$ IVDEP
        DO jc = 1,ict
          ic = idxt

          zw  =pw(jk)
          zg  =pgg(jk)

          zg3= 3._wp * zg
          zgamma1= (8._wp - zw * (5._wp + zg3)) * 0.25_wp
          zgamma2=  3._wp *(zw * (1._wp - zg )) * 0.25_wp
          zgamma3= (2._wp - zg3 * prmuz ) * 0.25_wp
          zzz=(1._wp - zg)**2
          zcrit = zw*zzz - zwcrit*(zzz - (1._wp - zw)*(zg **2))
        END DO
      ELSEIF (kmodts == 3) THEN
!IBM* ASSERT(NODEPS)
!CDIR NODEP,VOVERTAKE,VOB
        DO jc = 1,ict
          ic = idxt

          zw  =pw(jk)
          zg  =pgg(jk)

          zg3= 3._wp * zg
          zgamma1= zsr3 * (2._wp - zw * (1._wp + zg)) * 0.5_wp
          zgamma2= zsr3 * zw * (1._wp - zg ) * 0.5_wp
          zgamma3= (1._wp - zsr3 * zg * prmuz ) * 0.5_wp
          zzz=(1._wp - zg)**2
          zcrit = zw*zzz - zwcrit*(zzz - (1._wp - zw)*(zg **2))
        END DO
      ENDIF

      icc = 0
      icn = 0
!IBM* ASSERT(NODEPS)
      DO jc = 1,ict
        ic = idxt

        !-- RECOMPUTE ORIGINAL S.S.A. TO TEST FOR CONSERVATIVE SOLUTION
        !   ZTEMP=(1._wp - ZG)**2
        !   ZWO= ZW*ZTEMP/ (ZTEMP - (1._wp - ZW)*(ZG **2))

        !       ZWO= ZW / (1._wp - (1._wp - ZW) * (ZG / (1._wp - ZG))**2)
        !       IF (ZWO >= ZWCRIT) THEN

        IF (zcrit >= 0._wp) THEN
           icc = icc + 1
           idxc = ic
        ELSE
           icn = icn + 1
           idxn = ic
        END IF
      END DO

       !!-- conservative scattering

        !-- Homogeneous reflectance and transmittance

      IF (icc == icount) THEN ! use direct addressing
        DO ic = 1,icc
          zem2 = -MIN(ptau(jk) * prmuzi,500._wp)
        END DO
      ELSE
!IBM* ASSERT(NODEPS)
!CDIR NODEP,VOVERTAKE,VOB
        DO jc = 1,icc
          ic = idxc
          zem2 = -MIN(ptau(jk) * prmuzi,500._wp)
        END DO
      ENDIF

      zem2 = EXP(zem2)

      IF (icc == icount) THEN ! use direct addressing
!IBM* NOVECTOR
        DO ic = 1,icc
          za  = zgamma1 * prmuz
          za1 = za - zgamma3
          zgt = zgamma1 * ptau(jk)

          ! collimated beam

          ztemp=1.0_wp/(1._wp + zgt)
          pref(jk) = (zgt - za1 * (1._wp - zem2)) *ztemp
          ptra(jk) = 1._wp - pref(jk)

          ! isotropic incidence

          prefd(jk) = zgt *ztemp
          ptrad(jk) = 1._wp - prefd(jk)
        END DO
      ELSE
!IBM* NOVECTOR
!IBM* ASSERT(NODEPS)
!CDIR NODEP,VOVERTAKE,VOB
!DIR$ IVDEP
        DO jc = 1,icc
          ic = idxc

          za  = zgamma1 * prmuz
          za1 = za - zgamma3
          zgt = zgamma1 * ptau(jk)

          ! collimated beam

          ztemp=1.0_wp/(1._wp + zgt)
          pref(jk) = (zgt - za1 * (1._wp - zem2)) *ztemp
          ptra(jk) = 1._wp - pref(jk)

          ! isotropic incidence

          prefd(jk) = zgt *ztemp
          ptrad(jk) = 1._wp - prefd(jk)
        END DO
      ENDIF

     !-- non-conservative scattering

     !-- Homogeneous reflectance and transmittance
      IF (icn == icount) THEN ! use direct addressing
        DO ic = 1,icn
          zzz = zgamma1**2 - zgamma2**2
          zrk = SQRT(MAX(zzz,replog))

          zep1 = MIN(zrk * ptau(jk), 500._wp)
          zep2 = MIN(ptau(jk) * prmuzi,500._wp)
        END DO
      ELSE
!IBM* ASSERT(NODEPS)
!CDIR NODEP,VOVERTAKE,VOB
!DIR$ IVDEP
        DO jc = 1,icn
          ic = idxn
          zzz = zgamma1**2 - zgamma2**2
          zrk = SQRT(MAX(zzz,replog))

          zep1 = MIN(zrk * ptau(jk), 500._wp)
          zep2 = MIN(ptau(jk) * prmuzi,500._wp)
        END DO
      ENDIF

      zep1 = EXP(zep1)
      zep2 = EXP(zep2)
      zem1 = 1.0_wp/zep1
      zem2 = 1.0_wp/zep2

      IF (icn == icount) THEN ! use direct addressing
        DO ic = 1,icn
          zw  =pw(jk)

          zgamma4 = 1._wp - zgamma3

          za1 = zgamma1 * zgamma4 + zgamma2 * zgamma3
          za2 = zgamma1 * zgamma3 + zgamma2 * zgamma4

          zrp = zrk * prmuz
          zrp1 = 1._wp + zrp
          zrm1 = 1._wp - zrp
          zrk2 = 2._wp * zrk
          zrpp = 1._wp - zrp*zrp
          zrkg = zrk + zgamma1
          zr1  = zrm1 * (za2 + zrk * zgamma3)
          zr2  = zrp1 * (za2 - zrk * zgamma3)
          zr3  = zrk2 * (zgamma3 - za2 * prmuz )
          zr4  = zrpp * zrkg
          zr5  = zrpp * (zrk - zgamma1)
          zt1  = zrp1 * (za1 + zrk * zgamma4)
          zt2  = zrm1 * (za1 - zrk * zgamma4)
          zt3  = zrk2 * (zgamma4 + za1 * prmuz )
          ! zt4  = zr4
          ! zt5  = zr5
          ! GZ, 2015-06-08: another fix for potential division by zero
          zbeta = - zr5 / SIGN(MAX(zeps,ABS(zr4)),zr4)

          ! collimated beam

          ! GZ, 2014-07-03: provisional fix for potential division by zero
          zaux = zr4*zep1 + zr5*zem1
          zdenr = SIGN(MAX(zeps,ABS(zaux)),zaux)
          pref(jk) = zw  * (zr1*zep1 - zr2*zem1 - zr3*zem2) / zdenr

          ! zdent = zt4*zep1 + zt5*zem1
          zdent = zdenr
          ptra(jk) = zem2 * &
            & (1._wp - zw  * (zt1*zep1 - zt2*zem1 - zt3*zep2) / zdent)

          ! diffuse beam

          zemm = zem1*zem1
          zdend = 1._wp / ( (1._wp - zbeta*zemm ) * zrkg)
          prefd(jk) =  zgamma2 * (1._wp - zemm) * zdend
          ptrad(jk) =  zrk2*zem1*zdend

        END DO
      ELSE
!IBM* ASSERT(NODEPS)
!CDIR NODEP,VOVERTAKE,VOB
!DIR$ IVDEP
        DO jc = 1,icn
          ic = idxn
          zw  =pw(jk)

          zgamma4 = 1._wp - zgamma3

          za1 = zgamma1 * zgamma4 + zgamma2 * zgamma3
          za2 = zgamma1 * zgamma3 + zgamma2 * zgamma4

          zrp = zrk * prmuz
          zrp1 = 1._wp + zrp
          zrm1 = 1._wp - zrp
          zrk2 = 2._wp * zrk
          zrpp = 1._wp - zrp*zrp
          zrkg = zrk + zgamma1
          zr1  = zrm1 * (za2 + zrk * zgamma3)
          zr2  = zrp1 * (za2 - zrk * zgamma3)
          zr3  = zrk2 * (zgamma3 - za2 * prmuz )
          zr4  = zrpp * zrkg
          zr5  = zrpp * (zrk - zgamma1)
          zt1  = zrp1 * (za1 + zrk * zgamma4)
          zt2  = zrm1 * (za1 - zrk * zgamma4)
          zt3  = zrk2 * (zgamma4 + za1 * prmuz )
          ! zt4  = zr4
          ! zt5  = zr5
          ! GZ, 2015-06-08: another fix for potential division by zero
          zbeta = - zr5 / SIGN(MAX(zeps,ABS(zr4)),zr4)

          ! collimated beam

          ! GZ, 2014-07-03: provisional fix for potential division by zero
          zaux  = zr4*zep1 + zr5*zem1
          zdenr = SIGN(MAX(zeps,ABS(zaux)),zaux)
          pref(jk) = zw  * (zr1*zep1 - zr2*zem1 - zr3*zem2) / zdenr

          ! zdent  = zt4*zep1(jc) + zt5*zem1(jc)
          zdent = zdenr
          ptra(jk) = zem2 * &
            & (1._wp - zw  * (zt1*zep1 - zt2*zem1 - zt3*zep2) / zdent)

          ! diffuse beam

          zemm = zem1*zem1
          zdend = 1._wp / ( (1._wp - zbeta*zemm ) * zrkg)
          prefd(jk) =  zgamma2 * (1._wp - zemm) * zdend
          ptrad(jk) =  zrk2*zem1*zdend

        END DO
      ENDIF

      IF (icf == icount) THEN ! use direct addressing
        DO ic = 1,icf
          pref(jk) =0.0_wp
          ptra(jk) =1.0_wp
          prefd(jk)=0.0_wp
          ptrad(jk)=1.0_wp
        END DO
      ELSE IF (icf > 0) THEN
!IBM* ASSERT(NODEPS)
!CDIR NODEP,VOVERTAKE,VOB
        DO jc = 1,icf
          ic=idxf
          pref(jk) =0.0_wp
          ptra(jk) =1.0_wp
          prefd(jk)=0.0_wp
          ptrad(jk)=1.0_wp
        END DO
      ENDIF
    ENDDO

  END SUBROUTINE srtm_reftra

END MODULE mo_srtm
