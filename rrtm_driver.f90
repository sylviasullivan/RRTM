PROGRAM rrtm_driver
    USE rrtm
    IMPLICIT NONE
 
    ! >> sylvia_20200310
    LOGICAL  :: writeoutput
    INTEGER  ::         &
      &  day,                      & !< day associated with this profile, sylvia_20200305
      &  hour,                     & !< hour associated with this profile
      &  minute                      !< minute associated with this profile

    INTEGER,  PARAMETER :: pd =  12
    INTEGER,  PARAMETER :: rd = 307
    INTEGER,  PARAMETER :: dp = SELECTED_REAL_KIND(pd,rd) !< double precision
    INTEGER,  PARAMETER :: wp = dp
    INTEGER,  PARAMETER :: klev = 81 ! 90  >> sylvia_20200520
    REAL(wp), PARAMETER :: pi = 4 * ATAN(1.0_8) ! >> sylvia_20200602

    ! >> sylvia_20200520
    REAL, DIMENSION(9,81) :: ellingson
    INTEGER               :: i,indx,time
    REAL(wp)              :: qi_val
    REAL(wp)              :: sDeclin, delJ, J, L0, L1, L2, L3, lonn
    REAL(wp)              :: sidereal, LHA, taud
    ! << sylvia_20200520

    ! >> sylvia_20200310
    REAL(wp)  ::         &
      &  latt,                     & !< latitude associated with this profile, sylvia_20200305
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
      &  aotss(12),                & !< monthly climatology, sea salt AOT
      &  aotorg(12),               & !< monthly climatology, organic AOT
      &  aotbc(12),                & !< monthly climatology, black carbon AOT
      &  aotso4(12),               & !< monthly climatology, sulfate AOT
      &  aotdu(12)                   !< monthly climatology, dust AOT

    ! >> sylvia_20200310
    REAL(wp) ::                     &
      &  flx_lw_net(klev+1),        & !< net downward LW flux profile,
      &  flx_sw_net(klev+1),        & !< net downward SW flux profile,
      &  flx_lw_net_clr(klev+1),    & !< clrsky downward LW flux profile,
      &  flx_sw_net_clr(klev+1),    & !< clrsky downward SW flux profile,
      &  flx_uplw_sfc,              & !< sfc LW upward flux,
      &  flx_upsw_sfc,              & !< sfc SW upward flux,
      &  flx_uplw_sfc_clr,          & !< clrsky sfc LW upward flux,
      &  flx_upsw_sfc_clr             !< clrsky sfc SW upward flux,

    ! optional..
    REAL(wp) ::    &
      &  flx_upsw_toa                !< TOA SW upward flux,

    ! >> sylvia_20200406
    ! stackoverflow.com/questions/13843772/fortran-read-command-line-arguments-how-to-use-command-line-argument
    ! Pass in four inputs as alb_vis_dir, alb_nir_dir, alb_vis_dif, alb_nir_dif
    INTEGER  :: num_args, ix
    CHARACTER(LEN=12), DIMENSION(:), ALLOCATABLE :: args
    CHARACTER(LEN=10) :: file_tag

    ! >> sylvia_20200406
    !num_args = command_argument_count()
    !ALLOCATE(args(num_args))
    !DO ix = 1, num_args
    !    CALL get_command_argument(ix,args(ix))
    !END DO
    !READ(args(1),*) alb_vis_dir
    ! << sylvia_20200406

    ! >> sylvia_20200528
    num_args = command_argument_count()
    ALLOCATE(args(num_args))
    DO ix = 1, num_args
       CALL get_command_argument(ix,args(ix))
    END DO
    READ(args(1),*) file_tag
    READ(args(2),*) indx
    READ(args(3),*) qi_val
    !READ(args(3),*) time
    ! << sylvia_20200528

    ! Consider an ocean surface with low albedo at 7 pm in the tropics.
    alb_vis_dir = 0.05_wp
    alb_nir_dir = 0.05_wp
    alb_vis_dif = 0.05_wp
    alb_nir_dif = 0.05_wp
    day         = 8
    hour        = time
    minute      = 0
    latt        = 0.0872_wp  !  5 deg N to radians
    lonn        = 1.1344_wp  ! 65 deg E to radians
    zland       = 0._wp
    emis_rad    = 0.99_wp    ! Emissivities tend to be > 0.85 (sand has lowest)

    open(1,file='output/tropical_profile_ellingson_250m_formatted_top2bottom.txt')
    read(1,*) ellingson

    ! Full-level pressures in [Pa], Factor of 100 for hPa -> Pa
    pp_fl = ellingson(2,:)*100._wp

    ! Half-level pressures in [Pa]
    do i = 2, klev
       pp_hl(i) = (pp_fl(i-1) + pp_fl(i))/2.0_wp
    end do

    ! The first half-level is the first full-level minus the difference up to the second half-level.
    ! The last half-level is the last full-level plus the difference down to the second-to-last half-level.
    pp_hl(1)      = pp_fl(1) - (pp_hl(2) - pp_fl(1))
    pp_hl(klev+1) = pp_fl(klev) + (pp_fl(klev) - pp_hl(klev))

    ! Surface pressure in [Pa], Factor of 100 for hPa -> Pa
    pp_sfc = ellingson(2,1)*100._wp

    ! Full-level temperatures in [K]
    tk_fl = ellingson(3,:)

    ! Half-level temperatures in [K]
    do i = 2, klev
       tk_hl(i) = (tk_fl(i-1) + tk_fl(i))/2.0_wp
    end do
    tk_hl(1)      = tk_fl(1) - (tk_hl(2) - tk_fl(1))
    tk_hl(klev+1) = tk_fl(klev) + (tk_fl(klev) - tk_hl(klev))

    ! Surface temperature in [K]
    tk_sfc = ellingson(3,1)

    ! Specific humidity profile in [kg kg-1]
    xm_vap = ellingson(6,:)

    ! Specific cloud liquid profile in [kg kg-1]
    xm_liq(:)  = 0._wp

    ! Specific cloud ice water profile in [kg kg-1]
    xm_ice = (/ 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, &
              & 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, &
              & 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, &
              & 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, &
              & 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, &
              & 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, &
              & 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, &
              & 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp /)
    xm_ice(indx) = qi_val

    ! Cloud droplet number concentration [m-3], 0.1 cm-3 below
    ! For some reason I don't yet understand, we need a non-zero value below to avoid runtime errors.
    cdnc  = 0.1_wp

    ! Aerosol optical thicknesses
    ! Sea salt
    aotss = (/ 0.008853, 0.010085, 0.008805, 0.007506, 0.008948, 0.013439, 0.010637, &
             & 0.009285, 0.007854, 0.008478, 0.009137, 0.010712 /)
    ! Organic
    aotorg = (/ 0.02928 , 0.036032, 0.038364, 0.022213, 0.014781, 0.02086 , 0.040556, &
              & 0.037241, 0.028972, 0.024837, 0.019659, 0.026071 /)
    ! Black carbon
    aotbc = (/ 0.00266 , 0.003283, 0.003511, 0.00207 , 0.001312, 0.00208 , 0.004036, &
             & 0.003704, 0.002956, 0.002389, 0.001803, 0.002434 /)

    ! Sulfate
    aotso4 = (/ 0.023888, 0.020334, 0.020786, 0.019716, 0.025483, 0.019174, 0.028141, &
              & 0.03837 , 0.035887, 0.025932, 0.014476, 0.016929 /)

    ! Dust
    aotdu = (/ 0.00853149, 0.03173592, 0.0456236 , 0.03786197, 0.05928399, 0.2629021, &
             & 0.221159  , 0.2713769 , 0.13333197, 0.01799284, 0.01569963, 0.01770153 /)

    ! Take the cloud fraction to be binary for now and set it according to xm_ice + xm_liq above.
    ! It is assumed to be a dependent parameter.
    do i = 1, klev
       if (xm_liq(i) + xm_ice(i) > 1e-10) then
          cld_frc(i) = 1.0_wp
       else
          cld_frc(i) = 0.0_wp
       end if
    end do

    ! Solar zenith angle
    ! Solar declination as in http://naturalfrequency.com/Tregenza_Sharples/Daylight_Algorithms/algorithm_1_11.htm
    ! Month is set to August within the code. July 31st is the 212th day in a non-leap year.
    J        = 212 + day
    taud     = 2._wp*pi*(J - 1._wp)/360._wp
    sDeclin  = 0.006918_wp - 0.399912_wp*COS(taud) + 0.072057_wp*SIN(taud) - 0.006758_wp*COS(2*taud) + &
            &  0.000907_wp*SIN(2*taud) - 0.002697_wp*COS(3*taud) + 0.001480*SIN(3*taud)
    ! Local hour angle is the local sidereal time - right ascension.
    ! Local sidereal time as https://www.aa.quae.nl/en/reken/sterrentijd.html
    L0       = 99.967794687_wp
    L1       = 360.985647366_wp
    L2       = 2.9078e-13_wp
    L3       = -5.302e-22_wp
    ! Note that if the longitude changes, the factor of 6 below must also. It represents the time zone.
    ! I have also assumed 2017 is the year below (6210 days after 1 Jan 2000).
    delJ     = 6210._wp + J + (hour - 6._wp)/24._wp
    sidereal = L0 + L1*delJ +L2*delJ**2._wp + L3*delJ**3._wp - (-1._wp*lonn*180._wp/pi)
    ! Convert sidereal time to radians for use below.
    sidereal = sidereal*pi/180._wp
    LHA      = sidereal
    pmu0     = ACOS(SIN(latt)*SIN(sDeclin) + COS(latt)*COS(sDeclin)*COS(LHA))

    CALL rrtm_interface (               klev                             ,& 
    & day             ,hour            ,minute          ,latt            ,&
    & zland           ,pmu0            ,alb_vis_dir                      ,&
    & alb_nir_dir     ,alb_vis_dif     ,alb_nir_dif     ,emis_rad        ,&
    & pp_fl           ,pp_hl           ,pp_sfc                           ,&
    & tk_fl           ,tk_hl           ,tk_sfc                           ,&
    & xm_vap          ,xm_liq          ,xm_ice                           ,&
    & cdnc            ,cld_frc         ,aotss           ,aotorg          ,&
    & aotbc           ,aotso4          ,aotdu                            ,&
    ! >> sylvia_20200305
    ! Changing inputs zaeq1,zaeq2,zaeq3,zaeq4,zaeq5 to aot* variables.
    ! output
    & flx_lw_net      ,flx_sw_net      ,flx_lw_net_clr  ,flx_sw_net_clr  ,&
    & flx_uplw_sfc    ,flx_upsw_sfc    ,flx_uplw_sfc_clr,flx_upsw_sfc_clr,&
    ! optional output
    & flx_upsw_toa                                                        )

    ! Output the fluxes
    writeoutput = .true.
    if (writeoutput .eqv. .true.) then
        open(unit=1,file='lwflxatm-test' // trim(file_tag) // '.txt',status='new')
        write(1,*) flx_lw_net
        write(1,*) flx_lw_net_clr
        close(1)

        open(unit=2,file='swflxatm-test' // trim(file_tag) // '.txt',status='new')
        write(2,*) flx_sw_net
        write(2,*) flx_sw_net_clr
        close(2)

        open(unit=3,file='flxsfc-test' // trim(file_tag) // '.txt',status='new')
        write(3,*) flx_uplw_sfc
        write(3,*) flx_uplw_sfc_clr
        write(3,*) flx_upsw_sfc
        write(3,*) flx_upsw_sfc_clr
        close(3)
    endif

    !  Print the net LW, SW fluxes
    !write(*,*) 'flx_uplw_sfc',flx_uplw_sfc
    !write(*,*) 'flx_uplw_sfc_clr',flx_uplw_sfc_clr
    !write(*,*) 'flx_lw_net - flx_lw_net_clr',flx_lw_net - flx_lw_net_clr
    !write(*,*) 'flx_lw_net_clr',flx_lw_net_clr
    !write(*,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
    !write(*,*) 'flx_upsw_sfc',flx_upsw_sfc
    !write(*,*) 'flx_upsw_sfc_clr',flx_upsw_sfc_clr
    !write(*,*) 'flx_sw_net - flx_sw_net_clr',flx_sw_net - flx_sw_net_clr
    !write(*,*) 'flx_sw_net_clr',flx_sw_net_clr
    !write(*,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'

END PROGRAM rrtm_driver
